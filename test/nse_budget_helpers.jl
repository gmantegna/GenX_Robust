# Helpers shared by test_nse_budget.jl and nse_budget_stock_mes.jl (the subprocess that checks
# stock MacroEnergySolvers). `include` this file; it is not a test file.
#
# HiGHS threads: HiGHS >= 1.7 has ONE global scheduler per process, sized by the first solve.
# A later solve that asks for a different `threads` value gets "Option threads is set to 1 but
# global scheduler has already been initialized" -> OTHER_ERROR, and MES then spins in
# compute_conflict!. The Benders settings of the case use `threads: 1`, so EVERY HiGHS solve
# made by these tests uses threads = 1 (see `HIGHS_THREADS`, `solve_monolithic`).

using Test
using GenX
using JuMP, HiGHS
using Logging
using YAML, CSV, DataFrames

const MES = GenX.MacroEnergySolvers

# One thread setting for every HiGHS solve in the process; see the header.
const HIGHS_THREADS = 1
# If an earlier test file in the same process (runtests.jl) already sized the scheduler with
# the default thread count, reset it so that the threads = 1 solves below are accepted.
if isdefined(HiGHS, :Highs_resetGlobalScheduler)
    HiGHS.Highs_resetGlobalScheduler(1)
end

const CASE = joinpath(@__DIR__, "benders", "1_three_zones")
const SETTINGS = joinpath(CASE, "settings")

# $/MWh. Scaled down from the case's 50000 $/MWh; see header. Measured on this case (single
# segment, HiGHS): E* = 0 at VoLL = 200, 300, 500, 1000, 2000 $/MWh; E* = 263966 MWh/yr
# (0.23% of demand) at 100 $/MWh. The case is 24 modelled hours with cheap energy.
const VOLL_TEST = 100.0

# Same tolerance as test_benders_vs_monolithic.jl
const OBJECTIVE_RTOL = 1e-3
# Monolithic-vs-monolithic identities are LP-exact up to solver tolerance.
const LP_RTOL = 1e-6

# Known monolithic optimum of the UNMODIFIED case (from test_benders_vs_monolithic.jl).
const KNOWN_OPTIMUM_CASE1 = 4975.3803

quiet(f) = redirect_stdout(() -> with_logger(f, ConsoleLogger(stderr, Logging.Warn)), devnull)

function case_setup(overrides::Dict = Dict())
    setup = quiet() do
        GenX.configure_settings(joinpath(SETTINGS, "genx_settings.yml"),
            joinpath(SETTINGS, "output_settings.yml"))
    end
    setup["PrintModel"] = 0
    merge!(setup, overrides)
    return setup
end

scale(setup) = setup["ParameterScale"] == 1 ? GenX.ModelScalingFactor : 1.0

"Load the case; optionally collapse to one curtailment segment priced at `voll` (\$/MWh)."
function case_inputs(setup; single_segment::Bool = true, voll::Float64 = VOLL_TEST)
    inputs = quiet(() -> GenX.load_inputs(setup, CASE))
    if single_segment
        inputs["SEG"] = 1
        inputs["pC_D_Curtail"] = [voll / scale(setup)]
        inputs["pMax_D_Curtail"] = [1.0]
        inputs["Voll"] = [voll / scale(setup)]
    end
    return inputs
end

"Expected annual shed in MWh: sum_t omega[t] * sum_{s,z} vNSE[s,t,z]."
function expected_shed_mwh(EP, inputs, setup)
    nse = value.(EP[:vNSE])
    om = inputs["omega"]
    return scale(setup) *
           sum(om[t] * nse[s, t, z] for s in axes(nse, 1), t in axes(nse, 2), z in axes(nse, 3))
end

function solve_monolithic(setup, inputs)
    s = copy(setup)
    s["Benders"] = 0
    s["EnableJuMPStringNames"] = 1   # tests look variables up by name
    EP = quiet() do
        OPT = GenX.configure_solver(SETTINGS, HiGHS.Optimizer)
        EP = GenX.generate_model(s, inputs, OPT)
        set_silent(EP)
        set_attribute(EP, "threads", HIGHS_THREADS)
        optimize!(EP)
        EP
    end
    @assert termination_status(EP) == MOI.OPTIMAL
    return EP
end

"Logger that counts MES's feasibility-cut messages and swallows everything below Warn."
mutable struct FeasCutCounter <: AbstractLogger
    n::Int
end
Logging.min_enabled_level(::FeasCutCounter) = Logging.Info
Logging.shouldlog(::FeasCutCounter, args...) = true
function Logging.handle_message(l::FeasCutCounter, level, message, args...; kwargs...)
    if occursin("generating feasibility cut", string(message))
        l.n += 1
    end
    if level >= Logging.Warn
        println(stderr, "[$level] $message")
    end
end

function benders_setup(setup; convtol = 1e-4)
    s = copy(setup)
    s["Benders"] = 1
    sb = redirect_stdout(() -> GenX.configure_benders(joinpath(SETTINGS, "benders_settings.yml")), devnull)
    s = merge(s, sb)
    s["settings_path"] = SETTINGS
    s[:ConvTol] = convtol
    s[:ExpectFeasibleSubproblems] = false
    s[:Distributed] = false
    return s
end

function build_benders(s, inputs)
    redirect_stdout(devnull) do
        with_logger(ConsoleLogger(stderr, Logging.Warn)) do
            decomp = GenX.separate_inputs_subperiods(inputs)
            GenX.generate_benders_inputs(s, inputs, decomp, HiGHS.Optimizer)
        end
    end
end

"""
Run MES Benders on the in-memory case. Returns (results, n_feasibility_cuts,
initial budget allocation, expected shed in MWh at the incumbent, benders inputs).
"""
function solve_benders(setup, inputs; convtol = 1e-4, resolve_incumbent::Bool = true)
    s = benders_setup(setup; convtol)
    bi = build_benders(s, inputs)
    pp, subs, lvs = bi["planning_problem"], bi["subproblems"], bi["planning_variables_sub"]
    counter = FeasCutCounter(0)
    results = redirect_stdout(devnull) do
        with_logger(counter) do
            MES.benders(pp, subs, lvs, s)
        end
    end
    # Initial master iterate of the budget variables (first column of the history).
    names_all = name.(all_variables(pp))
    W = inputs["REP_PERIOD"]
    q0 = Float64[]
    if s["NSEBudget"] == 1
        for w in 1:W
            i = findfirst(==("vNSEbudget[$w]"), names_all)
            push!(q0, results.planning_sol_hist[i, 1])
        end
    end
    resolve_incumbent || return results, counter.n, q0, NaN, bi
    # Re-solve the subproblems at the incumbent (what run_genx_case_benders! does) and
    # measure the expected shed of the accepted solution.
    redirect_stdout(devnull) do
        with_logger(ConsoleLogger(stderr, Logging.Warn)) do
            MES.solve_subproblems(subs, results.planning_sol, true)
        end
    end
    shed = 0.0
    for sp in subs
        w = sp[:subproblem_index]
        nse = value.(sp[:model][:vNSE])
        Tw = ((w - 1) * inputs["hours_per_subperiod"] + 1):(w * inputs["hours_per_subperiod"])
        om = inputs["omega"][Tw]
        shed += sum(om[t] * nse[sg, t, z]
        for sg in axes(nse, 1), t in axes(nse, 2), z in axes(nse, 3))
    end
    return results, counter.n, q0, scale(setup) * shed, bi
end


budget_setup(target_mwh; remove_voll = 1) = case_setup(Dict("NSEBudget" => 1,
    "NSEBudgetTargetMWh" => target_mwh, "NSEBudgetRemoveVoLL" => remove_voll))

# Target of the end-to-end `run_genx_case!` test (unmodified case: 4 segments, VoLL 50000).
const E2E_TARGET_MWH = 100000.0

"""
Copy of the case in a temporary folder with `overrides` written to its genx_settings.yml, for
end-to-end `run_genx_case!` runs. The repo's case folder is never written to.
"""
function scratch_case(overrides::Dict)
    dir = joinpath(mktempdir(), "case")
    cp(CASE, dir)
    for d in filter(startswith("results"), readdir(dir))
        rm(joinpath(dir, d); recursive = true)
    end
    f = joinpath(dir, "settings", "genx_settings.yml")
    y = YAML.load_file(f)
    merge!(y, overrides)
    YAML.write_file(f, y)
    return dir
end

e2e_benders_case() = scratch_case(Dict("Benders" => 1, "NSEBudget" => 1,
    "NSEBudgetTargetMWh" => E2E_TARGET_MWH, "OverwriteResults" => 1, "PrintModel" => 0))
