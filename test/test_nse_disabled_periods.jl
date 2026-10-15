module TestNSEDisabledPeriods

# Tests of `NSEDisabledPeriods` (reimagining-ra #50, lane D): with the setting on, every hour of
# each representative period listed in `policies/NSE_disabled_periods.csv` gets a row
# `cNSEDisabled[s,t,z]: vNSE[s,t,z] == 0`, so the model must serve demand in those periods with
# capacity (built or existing) instead of shedding. Everything else is unchanged.
#
# Case: test/benders/1_three_zones (4 representative periods x 6 h, 3 zones, UNEQUAL Sub_Weights
# 2286/2556/2940/978, ParameterScale = 1), read unchanged from the repo and modified IN MEMORY
# after `load_inputs` (as test_nse_budget.jl does): one curtailment segment priced at a low VoLL
# so that the VoLL-form optimum sheds, and the demand of ONE period scaled up so that the shed is
# concentrated there ("demand exceeds capacity in one period only"). The loader itself is tested
# on a temporary copy of the case with the policy CSV written in.

include(joinpath(@__DIR__, "nse_budget_helpers.jl"))

# Period whose demand is scaled up and which the policy then disables: period 4, the one with
# the smallest weight (978 h), is the only period in which shedding at 100 $/MWh is cheaper than
# building (100 $/MWh × 978 h/yr per MW-yr of shed, below the cheapest capacity's annual cost;
# the other periods weigh 2286-2940 h, so demand there is met by building — measured in the lane
# D test job: with period 2 stressed the model built and did not shed). The mutant check of the
# lane report ("disable the wrong period") sets DISABLED_PERIOD = 2, a period that does not shed;
# the named test "disabling the shedding period removes its shed and builds capacity" then fails
# on `objective_value(EP_on) > Z_REF` (the row exists but binds nothing).
const STRESS_PERIOD = 4
const DISABLED_PERIOD = 4
const DEMAND_SCALE = 3.0
const VOLL_DISABLED_TEST = 100.0

"Per-period expected shed (MWh) and raw shed (MWh, unweighted) of a solved model."
function shed_by_period(EP, inputs, setup)
    nse = value.(EP[:vNSE])
    H = inputs["hours_per_subperiod"]
    W = inputs["REP_PERIOD"]
    raw = zeros(W)
    for w in 1:W, s in axes(nse, 1), t in ((w - 1) * H + 1):(w * H), z in axes(nse, 3)
        raw[w] += scale(setup) * nse[s, t, z]
    end
    return raw
end

function stressed_inputs(setup)
    inputs = case_inputs(setup; single_segment = true, voll = VOLL_DISABLED_TEST)
    H = inputs["hours_per_subperiod"]
    Tw = ((STRESS_PERIOD - 1) * H + 1):(STRESS_PERIOD * H)
    inputs["pD"] = copy(inputs["pD"])
    inputs["pD"][Tw, :] .*= DEMAND_SCALE
    return inputs
end

total_capacity_built(EP) = sum(value.(EP[:vCAP]))

const SETUP_REF = case_setup()
const INPUTS_REF = stressed_inputs(SETUP_REF)
const EP_REF = solve_monolithic(SETUP_REF, INPUTS_REF)
const Z_REF = objective_value(EP_REF)
const SHED_REF = shed_by_period(EP_REF, INPUTS_REF, SETUP_REF)
println("NSEDisabledPeriods tests: reference (VoLL=$(VOLL_DISABLED_TEST), period $STRESS_PERIOD demand x$DEMAND_SCALE): z = $Z_REF, shed by period (MWh) = $SHED_REF")

@testset "NSE disabled periods" begin
    @testset "T1 default off, setting validated" begin
        d = GenX.default_settings()
        @test d["NSEDisabledPeriods"] == 0
        @test !haskey(object_dictionary(EP_REF), :cNSEDisabled)
        # case_setup merges overrides after configure_settings has validated, so validate directly
        bad = merge(GenX.default_settings(), Dict{Any,Any}("NSEDisabledPeriods" => 2))
        @test_throws ErrorException GenX.validate_settings!(bad)
        ok = merge(GenX.default_settings(), Dict{Any,Any}("NSEDisabledPeriods" => 1))
        @test GenX.validate_settings!(ok) === nothing || true   # no throw
    end

    @testset "T2 loader reads policies/NSE_disabled_periods.csv" begin
        tmp = mktempdir()
        cp(CASE, joinpath(tmp, "case"); force = true)
        pol = joinpath(tmp, "case", "policies", "NSE_disabled_periods.csv")
        write(pol, "Rep_Period\n$(DISABLED_PERIOD)\n")
        s = case_setup(Dict("NSEDisabledPeriods" => 1))
        inp = quiet(() -> GenX.load_inputs(s, joinpath(tmp, "case")))
        @test inp["NSE_DISABLED_PERIODS"] == [DISABLED_PERIOD]
        H = inp["hours_per_subperiod"]
        @test GenX.nse_disabled_hours(inp) == collect(((DISABLED_PERIOD - 1) * H + 1):(DISABLED_PERIOD * H))
        # a Benders-style sliced input (SubPeriod set) maps to all/none of its own hours
        sub = copy(inp); sub["SubPeriod"] = DISABLED_PERIOD; sub["T"] = H
        @test GenX.nse_disabled_hours(sub) == collect(1:H)
        sub["SubPeriod"] = DISABLED_PERIOD + 1
        @test isempty(GenX.nse_disabled_hours(sub))
        # out-of-range index is an error
        write(pol, "Rep_Period\n$(inp["REP_PERIOD"] + 1)\n")
        @test_throws ErrorException quiet(() -> GenX.load_inputs(s, joinpath(tmp, "case")))
        rm(tmp; recursive = true, force = true)
    end

    @testset "reference sheds in the stressed period only" begin
        @test SHED_REF[STRESS_PERIOD] > 1.0
        for w in 1:length(SHED_REF)
            w == STRESS_PERIOD && continue
            @test SHED_REF[w] ≈ 0.0 atol = 1e-6
        end
    end

    @testset "disabling the shedding period removes its shed and builds capacity" begin
        s = case_setup(Dict("NSEDisabledPeriods" => 1))
        inp = copy(INPUTS_REF)
        inp["NSE_DISABLED_PERIODS"] = [DISABLED_PERIOD]
        EP_on = solve_monolithic(s, inp)
        shed_on = shed_by_period(EP_on, inp, s)
        println("  disabled period $DISABLED_PERIOD: z = $(objective_value(EP_on)), shed by period (MWh) = $shed_on, capacity built: ref $(total_capacity_built(EP_REF)) -> $(total_capacity_built(EP_on))")
        @test haskey(object_dictionary(EP_on), :cNSEDisabled)
        H = inp["hours_per_subperiod"]
        @test length(EP_on[:cNSEDisabled]) == inp["SEG"] * H * inp["Z"]
        @test shed_on[DISABLED_PERIOD] ≈ 0.0 atol = 1e-6
        # the model served the stressed demand with capacity instead of shed: it cost more, and
        # built more than the reference did
        @test objective_value(EP_on) > Z_REF * (1 + 1e-6)
        @test total_capacity_built(EP_on) > total_capacity_built(EP_REF) * (1 + 1e-6)
        # hours outside the listed period are untouched: no row there
        for t in 1:inp["T"]
            in_disabled = (DISABLED_PERIOD - 1) * H < t <= DISABLED_PERIOD * H
            @test (t in axes(EP_on[:cNSEDisabled], 2)) == in_disabled
        end
    end
end

end # module
