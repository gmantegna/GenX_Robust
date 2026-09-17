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

include(joinpath(@__DIR__, "nse_budget_helpers.jl"))

# ---------------------------------------------------------------------------
# Shared VoLL-form reference (single segment, VOLL_TEST): z_VoLL and E*
# ---------------------------------------------------------------------------
const SETUP_OFF = case_setup()
const INPUTS = case_inputs(SETUP_OFF)
const EP_VOLL = solve_monolithic(SETUP_OFF, INPUTS)
const Z_VOLL = objective_value(EP_VOLL)
const ESTAR_MWH = expected_shed_mwh(EP_VOLL, INPUTS, SETUP_OFF)
println("NSEBudget tests: VoLL form (VoLL=$(VOLL_TEST) \$/MWh, 1 segment): z_VoLL = $Z_VOLL, E* = $ESTAR_MWH MWh")

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
        # Weighted convention: omega[t] = Sub_Weights[w] / H for t in w, and the weights sum to
        # 8760, so sum_t omega[t]*NSE[t] (the capped quantity) is expected annual MWh.
        @test sum(INPUTS["Weights"]) ≈ 8760
        @test sum(INPUTS["omega"]) ≈ 8760
        @test length(unique(INPUTS["Weights"])) > 1      # unequal weights: a dropped omega shows
        for w in 1:INPUTS["REP_PERIOD"], t in 1:H
            @test INPUTS["omega"][(w - 1) * H + t] ≈ INPUTS["Weights"][w] / H
        end
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

    @testset "T3 planning cap: plain sum, all coefficients one" begin
        # vNSEbudget[w] already carries Sub_Weights[w] (weighted convention), so the cap must
        # NOT be weighted again: sum_w 1 * vNSEbudget[w] <= target / scale.
        target = 0.5 * ESTAR_MWH
        s = benders_setup(budget_setup(target))
        pp = build_benders(s, INPUTS)["planning_problem"]
        W = INPUTS["REP_PERIOD"]
        @test length(unique(INPUTS["Weights"])) > 1      # a pi_w-weighted cap would show
        cap = pp[:cNSEBudget_planning]
        @test length(pp[:vNSEbudget]) == W
        for w in 1:W
            @test normalized_coefficient(cap, pp[:vNSEbudget][w]) == 1.0
            @test lower_bound(pp[:vNSEbudget][w]) == 0.0
        end
        f = constraint_object(cap).func
        @test length(f.terms) == W                        # no other variable in the row
        @test all(==(1.0), values(f.terms))
        @test constraint_object(cap).set == MOI.LessThan(target / scale(s))
        @test normalized_rhs(cap) ≈ target / scale(s)
    end

    @testset "T5 settings validation" begin
        msg(f) =
            try
                f()
                "no error"
            catch e
                e isa ErrorException ? e.msg : "not an ErrorException: $(typeof(e))"
            end
        # target is mandatory, numeric and non-negative when the budget is on
        @test_throws ErrorException solve_monolithic(case_setup(Dict("NSEBudget" => 1)), INPUTS)
        @test occursin("NSEBudgetTargetMWh (expected annual unserved energy cap, MWh) is not set",
            msg(() -> solve_monolithic(case_setup(Dict("NSEBudget" => 1)), INPUTS)))
        @test_throws ErrorException solve_monolithic(budget_setup(-1.0), INPUTS)
        @test occursin("NSEBudgetTargetMWh must be a non-negative number (got -1.0)",
            msg(() -> solve_monolithic(budget_setup(-1.0), INPUTS)))
        @test occursin("NSEBudgetTargetMWh must be a non-negative number (got abc)",
            msg(() -> solve_monolithic(budget_setup("abc"), INPUTS)))
        # Benders: complete recourse does not hold, so ExpectFeasibleSubproblems: true is an error
        sb = benders_setup(budget_setup(0.5 * ESTAR_MWH))
        sb[:ExpectFeasibleSubproblems] = true
        @test_throws ErrorException build_benders(sb, INPUTS)
        @test occursin("NSEBudget = 1 requires ExpectFeasibleSubproblems: false",
            msg(() -> build_benders(sb, INPUTS)))
        # ... and without the budget the same setting is accepted
        sb_off = benders_setup(SETUP_OFF)
        sb_off[:ExpectFeasibleSubproblems] = true
        @test build_benders(sb_off, INPUTS) isa Dict
    end
end

@testset "MES patch: status test, gap check, guards" begin
    MOIU = MOI.Utilities
    # (termination status, primal status) -> does the subproblem count as solved?
    function solved(ts, ps)
        m = Model(() -> MOIU.MockOptimizer(MOIU.Model{Float64}()); add_bridges = false)
        @variable(m, x >= 0)
        MOIU.attach_optimizer(m)
        MOIU.set_mock_optimize!(unsafe_backend(m), mo -> MOIU.mock_optimize!(mo, ts, (ps, [0.0])))
        optimize!(m)
        @assert termination_status(m) == ts && primal_status(m) == ps
        return with_logger(() -> GenX.mes_has_feasible_point(m), NullLogger())
    end
    @test solved(MOI.OPTIMAL, MOI.FEASIBLE_POINT)
    @test !solved(MOI.INFEASIBLE, MOI.INFEASIBLE_POINT)          # the HiGHS warm-start case
    @test !solved(MOI.INFEASIBLE, MOI.NO_SOLUTION)
    @test !solved(MOI.INFEASIBLE_OR_UNBOUNDED, MOI.FEASIBLE_POINT)
    @test !solved(MOI.OTHER_ERROR, MOI.NO_SOLUTION)
    @test !solved(MOI.OPTIMAL, MOI.INFEASIBLE_POINT)
    @test solved(MOI.OPTIMAL, MOI.NEARLY_FEASIBLE_POINT)         # accepted, as stock MES does
    @test solved(MOI.ALMOST_OPTIMAL, MOI.NEARLY_FEASIBLE_POINT)
    @test solved(MOI.TIME_LIMIT, MOI.FEASIBLE_POINT)             # accepted, as stock MES does

    # negative gap beyond ConvTol is an error, within it a warning, positive gap is silent
    res(ub, lb) = (UB_hist = [Inf, ub], LB_hist = [0.0, lb], termination_status = "NEGATIVE GAP")
    @test_throws ErrorException GenX.check_benders_gap(res(5125.04, 6246.34), Dict(:ConvTol => 1e-3))
    @test_logs (:warn, r"negative gap within the tolerance") GenX.check_benders_gap(
        res(100.0, 100.0 + 1e-6), Dict(:ConvTol => 1e-3))
    @test_logs GenX.check_benders_gap(res(100.05, 100.0), Dict(:ConvTol => 1e-3))

    # load-time patch: in place, from GenX's file, replacing (not adding) the MES method
    @test get(ENV, "GENX_DISABLE_MES_PATCH", "") != "1"
    @test GenX.MES_FEASIBLE_POINT_SHIM_APPLIED[]
    @test pkgversion(MES) == GenX.MES_PATCHED_VERSION
    @test GenX.mes_source_has_defect()      # false once MES is fixed upstream: delete the patch
    meth = which(MES.solve_subproblem, (Model, NamedTuple, Vector{String}, Bool))
    @test basename(string(meth.file)) == "mes_feasible_point_shim.jl"
    @test length(methods(MES.solve_subproblem)) == 1
    @test GenX.mes_require_feasible_point!() === true    # idempotent public guard
end

# World age. MacroEnergySolvers 0.2.2 is patched by GenX at load time (GenX.__init__), so the
# patch is visible to every call made after `using GenX`, including calls made from inside a
# function. `guarded_benders` is the regression test of that: it calls the public guard and
# MES.benders within ONE function call, and it is the FIRST Benders solve of this process.
# (With the opt-in runtime shim of the first version of this PR the same call ran stock MES
# and ended with NEGATIVE GAP while the flag said the shim was applied.)
function guarded_benders(s, inputs)
    GenX.mes_require_feasible_point!()
    return solve_benders(s, inputs)
end

@testset "NSE budget, Benders (MES 0.2.2 patched at GenX load)" begin
    @testset "T3 Benders vs monolithic parity, budget mode" begin
        target = 0.5 * ESTAR_MWH
        s = budget_setup(target)
        EP = solve_monolithic(s, INPUTS)
        z_mono = objective_value(EP)
        @test GenX.MES_FEASIBLE_POINT_SHIM_APPLIED[]
        results, nfeas, q0, shed_bd, _ = guarded_benders(s, INPUTS)   # guard + solve in one function
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

    @testset "T6 run_genx_case! end to end (Benders, HiGHS, unmodified case)" begin
        # 4 segments, VoLL 50000 removed by the budget; target E2E_TARGET_MWH.
        s = budget_setup(E2E_TARGET_MWH)
        inputs_raw = case_inputs(s; single_segment = false)
        EP = with_logger(() -> solve_monolithic(s, inputs_raw), ConsoleLogger(stderr, Logging.Error))
        z_mono = objective_value(EP)
        @test z_mono≈4964.456 rtol=1e-6            # value measured by the reviewer of PR 9
        case = e2e_benders_case()
        redirect_stdout(devnull) do
            with_logger(ConsoleLogger(stderr, Logging.Error)) do
                GenX.run_genx_case!(case, HiGHS.Optimizer)
            end
        end
        out = joinpath(case, "results_benders")
        @test isdir(out)
        nse = CSV.read(joinpath(out, "nse.csv"), DataFrame; header = false)
        @test nse[3, 1] == "AnnualSum"
        nse_total = parse(Float64, string(nse[3, end]))
        nse_sum = sum(parse(Float64, string(x)) for x in nse[3, 2:(end - 1)])
        status = CSV.read(joinpath(out, "benders_convergence.csv"), DataFrame)
        ub = status[end, :UB]
        lb = status[end, :LB]
        println("T6 e2e: z_mono = $z_mono | results_benders: UB = $ub LB = $lb iters = $(nrow(status)) | nse.csv AnnualSum total = $nse_total MWh (sum of columns $nse_sum)")
        @test nse_total≈E2E_TARGET_MWH rtol=1e-5
        @test nse_sum≈nse_total rtol=1e-8
        @test all(parse(Float64, string(x)) >= -1e-6 for x in nse[3, 2:end])
        @test abs(ub - z_mono) / z_mono <= OBJECTIVE_RTOL
        @test lb <= z_mono * (1 + LP_RTOL)
    end

    @testset "T7 stock MES (subprocess, GENX_DISABLE_MES_PATCH=1)" begin
        # A fresh process with the escape hatch set, i.e. stock MacroEnergySolvers 0.2.2.
        script = joinpath(@__DIR__, "nse_budget_stock_mes.jl")
        cmd = addenv(ignorestatus(`$(Base.julia_cmd()) --project=$(Base.active_project()) $script`),
            "GENX_DISABLE_MES_PATCH" => "1")
        log = read(pipeline(cmd; stderr = devnull), String)
        kv = Dict(String(k) => String(v)
        for (k, v) in (split(l, "="; limit = 2) for l in split(log, '\n') if occursin("=", l)))
        for k in sort(collect(keys(kv)))
            println("T7 stock MES: $k = $(kv[k])")
        end
        @test occursin("STOCK DONE", log)
        @test get(kv, "PATCH_APPLIED", "") == "false"     # the escape hatch works
        # (a) The defect is still in MES. If these two become "unexpected pass", MES has been
        # fixed: delete src/benders/mes_feasible_point_shim.jl and this testset.
        @test_broken get(kv, "STOCK_STATUS", "") == "OPTIMAL"
        @test_broken parse(Float64, kv["STOCK_LB"]) <=
                     parse(Float64, kv["STOCK_ZMONO"]) * (1 + LP_RTOL)
        # (b) run_genx_case! must fail on a negative gap BEFORE writing any output
        @test startswith(get(kv, "RUNNER_ERROR", ""), "ErrorException")
        @test occursin("negative gap", lowercase(get(kv, "RUNNER_ERROR", "")))
        @test get(kv, "RUNNER_RESULTS_WRITTEN", "") == "false"
    end
end

end # module TestNSEBudget
