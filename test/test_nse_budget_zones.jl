module TestNSEBudgetZones

# `NSEBudgetZones` (reimagining-ra #173, lane panel-out): the expected-NSE budget restricted to a
# subset of zones. Case and in-memory modifications as test_nse_budget.jl (test/benders/1_three_zones,
# one curtailment segment at VOLL_TEST). With the setting:
#   * the budget row counts vNSE only in the listed zones (monolithic cNSEBudget and every
#     Benders subproblem's cNSEBudget);
#   * with NSEBudgetRemoveVoLL = 1 the NSE cost is removed in the listed zones ONLY; every other
#     zone keeps its segment prices (omega[t] * pC_D_Curtail[s] per MWh) in the objective.
# Absent = every zone: today's model, coefficient for coefficient (T1 compares the written LP).

include(joinpath(@__DIR__, "nse_budget_helpers.jl"))

const SETUP_OFF = case_setup()
const INPUTS = case_inputs(SETUP_OFF)
const Z = INPUTS["Z"]
const EP_VOLL = solve_monolithic(SETUP_OFF, INPUTS)

"Expected annual shed (MWh) of zone z."
zone_shed_mwh(EP, inputs, setup, z) = scale(setup) *
    sum(inputs["omega"][t] * value(EP[:vNSE][s, t, z]) for s in 1:inputs["SEG"], t in 1:inputs["T"])

const SHED_Z = [zone_shed_mwh(EP_VOLL, INPUTS, SETUP_OFF, z) for z in 1:Z]
const ZB = argmax(SHED_Z)                         # the budgeted zone: the one the VoLL form sheds most in
println("NSEBudgetZones tests: VoLL-form shed by zone (MWh) = $SHED_Z; budgeted zone $ZB")

zbudget_setup(target, zones; remove_voll = 1) = case_setup(Dict("NSEBudget" => 1,
    "NSEBudgetTargetMWh" => target, "NSEBudgetRemoveVoLL" => remove_voll, "NSEBudgetZones" => zones))

function build_monolithic(setup, inputs)
    s = copy(setup)
    s["Benders"] = 0
    s["EnableJuMPStringNames"] = 1
    return quiet() do
        GenX.generate_model(s, inputs, GenX.configure_solver(SETTINGS, HiGHS.Optimizer))
    end
end

lp_text(EP) = (f = tempname() * ".lp"; write_to_file(EP, f); read(f, String))

@testset "NSE budget zones" begin
    @testset "Z1 absent = all zones listed = today's model" begin
        target = 0.5 * sum(SHED_Z)
        EPa = build_monolithic(budget_setup(target), INPUTS)
        EPb = build_monolithic(zbudget_setup(target, collect(1:Z)), INPUTS)
        EPc = build_monolithic(zbudget_setup(target, reverse(collect(1:Z))), INPUTS)
        @test lp_text(EPa) == lp_text(EPb)
        @test lp_text(EPa) == lp_text(EPc)
        for EP in (EPa, EPb)
            set_silent(EP); set_attribute(EP, "threads", HIGHS_THREADS); optimize!(EP)
        end
        @test objective_value(EPa) == objective_value(EPb)
    end

    @testset "Z2 monolithic: row over the listed zones, VoLL elsewhere" begin
        @test SHED_Z[ZB] > 1.0
        target = 0.5 * SHED_Z[ZB]
        s = zbudget_setup(target, [ZB])
        EP = build_monolithic(s, INPUTS)
        row = EP[:cNSEBudget]
        obj = objective_function(EP)
        for t in (1, INPUTS["T"]), z in 1:Z
            v = EP[:vNSE][1, t, z]
            if z == ZB
                @test normalized_coefficient(row, v) ≈ INPUTS["omega"][t]
                @test coefficient(obj, v) == 0.0
            else
                @test normalized_coefficient(row, v) == 0.0
                @test coefficient(obj, v) ≈ s["ObjScale"] * INPUTS["omega"][t] * INPUTS["pC_D_Curtail"][1]
            end
        end
        set_silent(EP); set_attribute(EP, "threads", HIGHS_THREADS); optimize!(EP)
        @test termination_status(EP) == MOI.OPTIMAL
        shed_b = zone_shed_mwh(EP, INPUTS, s, ZB)
        println("Z2: target $target MWh on zone $ZB; shed there $shed_b MWh; z = $(objective_value(EP))")
        @test shed_b <= target * (1 + LP_RTOL)
        # the system-wide budget at the same target is a different model (it also counts the other zones)
        EPs = solve_monolithic(budget_setup(target), INPUTS)
        @test !isapprox(objective_value(EP), objective_value(EPs); rtol = 1e-6)
    end

    @testset "Z3 Benders: subproblem rows over the listed zones; Benders = monolith" begin
        target = 0.5 * SHED_Z[ZB]
        s = zbudget_setup(target, [ZB])
        bi = build_benders(benders_setup(s), INPUTS)
        H = INPUTS["hours_per_subperiod"]
        for sp in bi["subproblems"]
            w, m = sp[:subproblem_index], sp[:model]
            for t in (1, H), z in 1:Z
                c = normalized_coefficient(m[:cNSEBudget], m[:vNSE][1, t, z])
                z == ZB ? (@test c ≈ INPUTS["omega"][(w - 1) * H + t]) : (@test c == 0.0)
            end
        end
        z_mono = objective_value(solve_monolithic(s, INPUTS))
        res, _, _, _, _ = solve_benders(s, INPUTS; resolve_incumbent = false)
        println("Z3: monolith $z_mono, Benders UB $(res.UB_hist[end]) LB $(res.LB_hist[end])")
        @test res.UB_hist[end] ≈ z_mono rtol = OBJECTIVE_RTOL
        @test res.LB_hist[end] <= z_mono * (1 + 1e-6)
    end

    @testset "Z4 validation" begin
        t = 0.5 * SHED_Z[ZB]
        for bad in (Int[], [0], [Z + 1], [1, 1], "1", [1.5])
            @test_throws ErrorException build_monolithic(zbudget_setup(t, bad), INPUTS)
        end
        # NSEBudget: 0 makes the setting inert
        s0 = case_setup(Dict("NSEBudget" => 0, "NSEBudgetZones" => [ZB]))
        @test objective_value(solve_monolithic(s0, INPUTS)) == objective_value(EP_VOLL)
    end
end

end # module
