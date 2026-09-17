module TestNSEBudget

# Tests of the hard expected-unserved-energy (EUE) budget policy (`NSEBudget`),
# monolithic and Benders. Case: test/benders/1_three_zones (4 representative periods x 6 h,
# 3 zones, UNEQUAL Sub_Weights 2286/2556/2940/978 summing to 8760, ParameterScale = 1).
#
# The case is read from the repo unchanged; two modifications are made IN MEMORY after
# `load_inputs` (nothing is written to the case folder):
#   * the four demand-curtailment segments are collapsed to ONE segment (full demand
#     curtailable), because the VoLL/budget identity z_VoLL = z_budget + VoLL*E holds for a
#     single price of shed energy, not for a segment ladder under one cap;
#   * VoLL is scaled down from 50000 $/MWh to `VOLL_TEST` so that the VoLL-form optimum
#     sheds a strictly positive expected energy.
#
# Units: with ParameterScale = 1 the model works in GW / GWh / M$. All quantities named
# `*_mwh` below are converted to MWh; objectives are left in model units (M$).

using Test
using GenX
using JuMP, HiGHS
using Logging

const MES = GenX.MacroEnergySolvers

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
function solve_benders(setup, inputs; convtol = 1e-4)
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

# ---------------------------------------------------------------------------
# Shared VoLL-form reference (single segment, VOLL_TEST): z_VoLL and E*
# ---------------------------------------------------------------------------
const SETUP_OFF = case_setup()
const INPUTS = case_inputs(SETUP_OFF)
const EP_VOLL = solve_monolithic(SETUP_OFF, INPUTS)
const Z_VOLL = objective_value(EP_VOLL)
const ESTAR_MWH = expected_shed_mwh(EP_VOLL, INPUTS, SETUP_OFF)
println("NSEBudget tests: VoLL form (VoLL=$(VOLL_TEST) \$/MWh, 1 segment): z_VoLL = $Z_VOLL, E* = $ESTAR_MWH MWh")

budget_setup(target_mwh; remove_voll = 1) = case_setup(Dict("NSEBudget" => 1,
    "NSEBudgetTargetMWh" => target_mwh, "NSEBudgetRemoveVoLL" => remove_voll))

@testset "NSE budget" begin
    @testset "T1 budget off: nothing changes" begin
        # (a) defaults exist and are off
        d = GenX.default_settings()
        @test d["NSEBudget"] == 0
        @test d["NSEBudgetRemoveVoLL"] == 1
        # (b) the unmodified case (4 segments, VoLL 50000) still reproduces the known optimum
        inputs_raw = case_inputs(SETUP_OFF; single_segment = false)
        EP_raw = solve_monolithic(SETUP_OFF, inputs_raw)
        @test objective_value(EP_raw)≈KNOWN_OPTIMUM_CASE1 rtol=LP_RTOL
        @test !haskey(object_dictionary(EP_raw), :cNSEBudget)
        @test !haskey(object_dictionary(EP_raw), :vNSEbudget)
        # (c) NSEBudget: 0 with the other keys set is inert: same model size, same optimum
        s0 = case_setup(Dict("NSEBudget" => 0, "NSEBudgetTargetMWh" => 1.0,
            "NSEBudgetRemoveVoLL" => 1))
        EP0 = solve_monolithic(s0, inputs_raw)
        @test num_variables(EP0) == num_variables(EP_raw)
        @test num_constraints(EP0; count_variable_in_set_constraints = true) ==
              num_constraints(EP_raw; count_variable_in_set_constraints = true)
        @test objective_value(EP0) == objective_value(EP_raw)
        # (d) Benders master with the budget off has no budget variable
        bi = build_benders(benders_setup(SETUP_OFF), inputs_raw)
        @test !any(startswith("vNSEbudget"), name.(all_variables(bi["planning_problem"])))
    end

    @testset "T2 monolithic cap binds" begin
        @test ESTAR_MWH > 1.0   # the VoLL-form reference sheds
        target = 0.5 * ESTAR_MWH
        for remove_voll in (1, 0)
            s = budget_setup(target; remove_voll)
            EP = solve_monolithic(s, INPUTS)
            shed = expected_shed_mwh(EP, INPUTS, s)
            println("T2 remove_voll=$remove_voll: target = $target MWh, shed = $shed MWh, z = $(objective_value(EP))")
            @test shed <= target * (1 + LP_RTOL)
            @test shed >= target * (1 - LP_RTOL)   # binding: target is below E*
            if remove_voll == 0
                # priced variant at a binding cap costs more than the VoLL form
                @test objective_value(EP) > Z_VOLL
            end
        end
        # a slack cap (2 E*) in the priced variant reproduces the VoLL form
        s = budget_setup(2 * ESTAR_MWH; remove_voll = 0)
        EP = solve_monolithic(s, INPUTS)
        @test objective_value(EP)≈Z_VOLL rtol=LP_RTOL
    end

    @testset "T3 subproblem budget row carries omega" begin
        s = benders_setup(budget_setup(0.5 * ESTAR_MWH))
        bi = build_benders(s, INPUTS)
        H = INPUTS["hours_per_subperiod"]
        for sp in bi["subproblems"]
            w = sp[:subproblem_index]
            m = sp[:model]
            @test "vNSEbudget[$w]" in sp[:linking_variables_sub]
            row = m[:cNSEBudget]
            for t in (1, H), z in 1:INPUTS["Z"]
                @test normalized_coefficient(row, m[:vNSE][1, t, z]) ≈
                      INPUTS["omega"][(w - 1) * H + t]
            end
            @test normalized_coefficient(row, variable_by_name(m, "vNSEbudget[$w]")) ≈ -1.0
        end
    end

    @testset "T3 stock MES 0.2.2 + HiGHS: infeasible subproblem taken for solved" begin
        # Documents the MES defect that GenX.mes_require_feasible_point!() works around (see its
        # docstring). Must run BEFORE the shim is applied. If MES is fixed these turn into
        # "unexpected pass" errors: then delete the shim and this testset.
        target = 0.5 * ESTAR_MWH
        s = budget_setup(target)
        z_mono = objective_value(solve_monolithic(s, INPUTS))
        results, nfeas, _, _, _ = solve_benders(s, INPUTS)
        println("T3 stock MES: z_mono = $z_mono | UB = $(results.UB_hist[end]) LB = $(results.LB_hist[end]) iters = $(length(results.UB_hist)) status = $(results.termination_status) | feasibility cuts = $nfeas")
        @test_broken results.termination_status == "OPTIMAL"
        @test_broken results.LB_hist[end] <= z_mono * (1 + LP_RTOL)
    end
end

# World age: the shim redefines a MES method, so it has to be applied in its own top-level
# statement, between testsets, for the redefinition to be visible to the tests below.
GenX.mes_require_feasible_point!()

@testset "NSE budget, Benders (MES feasible-point shim applied)" begin
    @testset "T3 Benders vs monolithic parity, budget mode" begin
        target = 0.5 * ESTAR_MWH
        s = budget_setup(target)
        EP = solve_monolithic(s, INPUTS)
        z_mono = objective_value(EP)
        results, nfeas, q0, shed_bd, _ = solve_benders(s, INPUTS)
        ub, lb = results.UB_hist[end], results.LB_hist[end]
        println("T3: z_mono = $z_mono | Benders UB = $ub LB = $lb iters = $(length(results.UB_hist)) status = $(results.termination_status) | feasibility cuts = $nfeas | q0 = $q0 | incumbent shed = $shed_bd MWh (target $target)")
        @test results.termination_status == "OPTIMAL"
        @test abs(ub - lb) / max(abs(ub), 1.0) <= OBJECTIVE_RTOL
        @test abs(ub - z_mono) / max(abs(z_mono), 1.0) <= OBJECTIVE_RTOL
        @test lb <= z_mono * (1 + LP_RTOL)               # LB is a lower bound
        @test all(iszero, q0)                            # run starts from a zero allocation
        @test nfeas >= 1                                 # ... and needs feasibility cuts
        @test shed_bd <= target * (1 + 1e-5)             # accepted UB satisfies the cap
    end

    @testset "T3 feasibility cut: excludes the trial, keeps the monolithic optimum" begin
        target = 0.5 * ESTAR_MWH
        s = budget_setup(target)
        EP = solve_monolithic(s, INPUTS)
        sb = benders_setup(s)
        bi = build_benders(sb, INPUTS)
        pp, subs, lvs = bi["planning_problem"], bi["subproblems"], bi["planning_variables_sub"]
        W = INPUTS["REP_PERIOD"]
        H = INPUTS["hours_per_subperiod"]
        sol = redirect_stdout(devnull) do
            with_logger(ConsoleLogger(stderr, Logging.Error)) do
                MES.add_slacks_to_subproblems!(subs)
                MES.add_approximate_variable_cost!(pp, Int(W))
                # trial point: master optimum with no cuts, budget allocation forced to zero
                for w in 1:W
                    fix(pp[:vNSEbudget][w], 0.0; force = true)
                end
                planning_sol, _ = MES.solve_planning_problem(pp, name.(all_variables(pp)))
                for w in 1:W
                    unfix(pp[:vNSEbudget][w])
                    set_lower_bound(pp[:vNSEbudget][w], 0.0)
                end
                subop = MES.solve_subproblems(subs, planning_sol, false)
                MES.update_planning_problem_multi_cuts!(pp, subop, planning_sol, lvs)
                (planning_sol, subop)
            end
        end
        planning_sol, subop = sol
        infeasible_w = [w for w in keys(subop) if subop[w].theta_coeff == 0]
        println("T3 feasibility: infeasible subproblems at the zero-budget trial = $(sort(infeasible_w)), slack* = $([subop[w].op_cost for w in sort(infeasible_w)])")
        @test !isempty(infeasible_w)
        @test all(subop[w].op_cost > 1e-6 for w in infeasible_w)
        # (i) the trial point is cut off: fixing every planning variable to it is infeasible
        trial = Dict(v => planning_sol.values[name(v)]
        for v in all_variables(pp) if !startswith(name(v), "vTHETA") && name(v) != "vZERO")
        for (v, val) in trial
            fix(v, val; force = true)
        end
        set_silent(pp)
        optimize!(pp)
        @test termination_status(pp) in (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
        # (ii) the monolithic optimum (x*, q*) survives the cuts
        nse = value.(EP[:vNSE])
        qstar = [sum(INPUTS["omega"][t] * nse[sg, t, z] for sg in axes(nse, 1),
                 t in ((w - 1) * H + 1):(w * H), z in axes(nse, 3)) for w in 1:W]
        n_matched = 0
        for v in keys(trial)
            nm = name(v)
            if startswith(nm, "vNSEbudget")
                w = parse(Int, nm[(length("vNSEbudget[") + 1):(end - 1)])
                fix(v, qstar[w]; force = true)
                n_matched += 1
            else
                vm = variable_by_name(EP, nm)
                @assert !isnothing(vm) "planning variable $nm not found in the monolithic model"
                fix(v, value(vm); force = true)
                n_matched += 1
            end
        end
        @test n_matched == length(trial)
        optimize!(pp)
        @test termination_status(pp) == MOI.OPTIMAL
    end

    @testset "T4 VoLL/budget identity" begin
        s = budget_setup(ESTAR_MWH)
        EP = solve_monolithic(s, INPUTS)
        z_budget = objective_value(EP)
        shed = expected_shed_mwh(EP, INPUTS, s)
        # VoLL * E* in model units (M\$ when ParameterScale = 1: \$/MWh * MWh / 1e6)
        offset = VOLL_TEST * ESTAR_MWH / scale(s)^2
        println("T4 monolithic: z_VoLL = $Z_VOLL, z_budget = $z_budget, VoLL*E* = $offset, z_budget + VoLL*E* = $(z_budget + offset), shed = $shed MWh, E* = $ESTAR_MWH MWh")
        @test Z_VOLL≈z_budget + offset rtol=LP_RTOL
        @test shed <= ESTAR_MWH * (1 + LP_RTOL)
        # the same identity through the Benders path
        results, nfeas, _, shed_bd, _ = solve_benders(s, INPUTS)
        ub = results.UB_hist[end]
        println("T4 Benders: UB = $ub, UB + VoLL*E* = $(ub + offset), feasibility cuts = $nfeas, incumbent shed = $shed_bd MWh")
        @test abs(ub + offset - Z_VOLL) / Z_VOLL <= OBJECTIVE_RTOL
        @test shed_bd <= ESTAR_MWH * (1 + 1e-5)
    end

    @testset "T5 settings validation" begin
        # target is mandatory and non-negative when the budget is on
        @test_throws Exception solve_monolithic(case_setup(Dict("NSEBudget" => 1)), INPUTS)
        @test_throws Exception solve_monolithic(budget_setup(-1.0), INPUTS)
    end
end

end # module TestNSEBudget
