module TestOpsFloorsBudgets

# Hourly Pmin, energy budgets, their slacks, the budgeted-hydro reservoir gate and the hourly hydro
# maximum (reimagining-ra #204, spec mwe/specs/base_ops_port.md, PR item 2). Hand values from the
# spec's section 9 (omega = 1, VoLL 1000 $/MWh, G = generator at cost 0, T = thermal at 10 $/MWh).

include(joinpath(@__DIR__, "ops_rows_helpers.jl"))

const ATOL = 1e-6

gen1h(; kw...) = ops_case(; T = 1, demand = [50.0], vre = [vre_row("G", 1, 100)],
    thermal = [thermal_row("TH", 1, 100, 10)], kw...)

# T3 family: 24 h, hydro 100 MW with budget column e, thermal T 100 MW at 10 $/MWh (optional)
function hydro24(; demand, e, therm = 100, var_h = 1.0, pmin = nothing, extra_hydro = String[],
        extra_var = Dict(), kw...)
    hydro = [hydro_row("Hy", 1, 100); extra_hydro]
    thermal = therm > 0 ? [thermal_row("TH", 1, therm, 10)] : String[]
    var = merge(Dict("Hy" => var_h isa Vector ? var_h : fill(var_h, length(demand))), extra_var)
    return ops_case(; T = length(demand), demand = demand, hydro = hydro, thermal = thermal,
        budget = Dict("Hy" => e isa Vector ? e : fill(e, length(demand))), variability = var,
        pmin = pmin, kw...)
end

@testset "Operational rows: floors, budgets, hydro maximum" begin
    @testset "T1 Pmin: T >= 30, objective 300" begin
        dir = gen1h(pmin = Dict("TH" => [0.3]))
        EP, inputs = ops_solve(dir; overrides = Dict("Hourly_Pmin" => 1))
        @test termination_status(EP) == MOI.OPTIMAL
        @test objective_value(EP)≈300 atol=ATOL
        @test vp(EP, inputs, "TH")[1]≈30 atol=ATOL
    end

    @testset "T2 floor slack: u = 10, objective 100,200; hard row infeasible" begin
        dir = gen1h(pmin = Dict("TH" => [0.3]), variability = Dict("TH" => [0.2]))
        EP, inputs = ops_solve(dir;
            overrides = Dict("Hourly_Pmin" => 1, "OperationalFloorPenalty" => 10_000))
        @test termination_status(EP) == MOI.OPTIMAL
        @test objective_value(EP)≈100_200 atol=ATOL
        @test sum(value.(EP[:vFloorSlack]))≈10 atol=ATOL
        @test value(EP[:eTotalCFloorSlack])≈100_000 atol=ATOL
        EPh, _ = ops_solve(dir; overrides = Dict("Hourly_Pmin" => 1))
        @test termination_status(EPh) in (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
    end

    @testset "T2b tie: NSE 270 at P = 10,000 and 1,001; P = 1,000 refused" begin
        dem = [fill(30.0, 23); 300.0]
        mk() = ops_case(; T = 24, demand = dem, hydro = [hydro_row("Hy", 1, 100)],
            vre = [vre_row("G", 1, 100)], budget = Dict("Hy" => fill(0.3, 24)),
            pmin = Dict("Hy" => fill(0.3, 24)),
            variability = Dict("G" => [fill(1.0, 23); 0.0], "Hy" => fill(1.0, 24)))
        dir = mk()
        EPh, _ = ops_solve(dir; overrides = Dict("Hourly_Pmin" => 1))
        @test nse_mwh(EPh)≈270 atol=ATOL
        for P in (10_000, 1_001)
            EP, _ = ops_solve(dir; overrides = Dict("Hourly_Pmin" => 1, "OperationalFloorPenalty" => P))
            @test termination_status(EP) == MOI.OPTIMAL
            @test nse_mwh(EP)≈270 atol=ATOL
            @test sum(value.(EP[:vFloorSlack]))≈0 atol=ATOL
        end
        @test occursin("must be a number above VoLL",
            errmsg_ops(() -> ops_solve(dir;
                overrides = Dict("Hourly_Pmin" => 1, "OperationalFloorPenalty" => 1_000))))
    end

    @testset "T3 budget: hydro 1,200 MWh, T 240 MWh, objective 2,400" begin
        dir = hydro24(demand = fill(60.0, 24), e = 0.5)
        EP, inputs = ops_solve(dir)
        @test termination_status(EP) == MOI.OPTIMAL
        @test objective_value(EP)≈2_400 atol=ATOL
        @test sum(vp(EP, inputs, "Hy"))≈1_200 atol=ATOL
        @test sum(vp(EP, inputs, "TH"))≈240 atol=ATOL
    end

    @testset "T3m hydro maximum: key 0 objective 6,000; key 1 NSE 60, objective 65,400" begin
        dem = [fill(100.0, 6); fill(40.0, 18)]
        dir = hydro24(demand = dem, e = 0.3, therm = 50, var_h = 0.4)
        EP0, inputs = ops_solve(dir; overrides = Dict("HydroBudgetMaxProfile" => 0))
        @test nse_mwh(EP0)≈0 atol=ATOL
        @test sum(vp(EP0, inputs, "TH"))≈600 atol=ATOL
        @test objective_value(EP0)≈6_000 atol=ATOL
        EP1, _ = ops_solve(dir; overrides = Dict("HydroBudgetMaxProfile" => 1))
        @test nse_mwh(EP1)≈60 atol=ATOL
        @test sum(vp(EP1, inputs, "TH"))≈540 atol=ATOL
        @test objective_value(EP1)≈65_400 atol=ATOL
        @test maximum(vp(EP1, inputs, "Hy")) <= 40 + ATOL
    end

    @testset "T3w window: 24 h windows 3,120,000; one 48 h window 480,000" begin
        dir = hydro24(demand = fill(50.0, 48), e = [fill(0.2, 24); fill(0.6, 24)], therm = 0)
        EP24, _ = ops_solve(dir;
            overrides = Dict("EnergyBudgetWindowHours" => 24, "OperationalFloorPenalty" => 10_000))
        @test termination_status(EP24) == MOI.OPTIMAL
        @test nse_mwh(EP24)≈720 atol=ATOL
        @test sum(value.(EP24[:vBudgetSpill]))≈240 atol=ATOL
        @test objective_value(EP24)≈3_120_000 atol=1e-4   # a window offset by 1 h: 2,680,000
        EP48, _ = ops_solve(dir; overrides = Dict("OperationalFloorPenalty" => 10_000))
        @test nse_mwh(EP48)≈480 atol=ATOL
        @test objective_value(EP48)≈480_000 atol=1e-4
    end

    @testset "T-spill: v = 720; variant (no upward slack) NSE 56" begin
        dir = hydro24(demand = fill(20.0, 24), e = 0.5, therm = 0)
        EP, inputs = ops_solve(dir; overrides = Dict("OperationalFloorPenalty" => 10_000))
        @test sum(vp(EP, inputs, "Hy"))≈480 atol=ATOL
        @test sum(value.(EP[:vBudgetSpill]))≈720 atol=ATOL
        @test value(EP[:eTotalCBudgetSpill])≈720 * 10_000 atol=1e-4
        @test objective_value(EP)≈7_200_000 atol=1e-4
        dirv = hydro24(demand = [fill(0.0, 23); 80.0], e = 0.01, therm = 0)
        EPv, inputsv = ops_solve(dirv; overrides = Dict("OperationalFloorPenalty" => 10_000))
        @test sum(vp(EPv, inputsv, "Hy"))≈24 atol=ATOL
        @test nse_mwh(EPv)≈56 atol=ATOL
        @test objective_value(EPv)≈56_000 atol=1e-4
    end

    @testset "T4 hydro without a budget column keeps the reservoir rows" begin
        dir = hydro24(demand = fill(60.0, 24), e = 0.5, extra_hydro = [hydro_row("Hy2", 1, 50)],
            extra_var = Dict("Hy2" => fill(0.1, 24)))
        EP, inputs = ops_solve(dir)
        g1 = findfirst(==("Hy"), inputs["RESOURCE_NAMES"])
        g2 = findfirst(==("Hy2"), inputs["RESOURCE_NAMES"])
        @test g2 in axes(EP[:cHydroReservoirInterior], 1)
        @test !(g1 in axes(EP[:cHydroReservoirInterior], 1))
        @test g2 in axes(EP[:cHydroReservoirStart], 1)
        @test !(g1 in axes(EP[:cHydroReservoirStart], 1))
        # Hy2 inflow 5 MWh/h = 120 MWh, Hy budget 1,200 MWh, demand 1,440: T 120 MWh
        @test objective_value(EP)≈1_200 atol=ATOL
    end

    @testset "T11 Benders subproblem rows of period 2 read the sliced floor and budget" begin
        prof = [fill(0.3, 24); fill(0.6, 24)]
        dir = ops_case(; T = 48, periods = 2, demand = fill(50.0, 48),
            thermal = [thermal_row("TH", 1, 100, 10)], hydro = [hydro_row("Hy", 1, 100)],
            pmin = Dict("TH" => prof), budget = Dict("Hy" => prof))
        setup = ops_setup(dir; overrides = Dict("Hourly_Pmin" => 1, "Benders" => 1))
        inputs = ops_inputs(dir, setup)
        gT = findfirst(==("TH"), inputs["RESOURCE_NAMES"])
        gH = findfirst(==("Hy"), inputs["RESOURCE_NAMES"])
        parts = GenX.separate_inputs_subperiods(inputs)
        OPT = optimizer_with_attributes(HiGHS.Optimizer)
        EP2 = quiet_ops(() -> GenX.generate_operation_subproblem(setup, parts[2], OPT))
        @test all(normalized_rhs(EP2[:cHourlyPmin][gT, t])≈60 for t in 1:24)     # unsliced: 30
        @test normalized_rhs(EP2[:cHourlyBudget][1, gH])≈1_440                       # unsliced: 720
        EP1 = quiet_ops(() -> GenX.generate_operation_subproblem(setup, parts[1], OPT))
        @test normalized_rhs(EP1[:cHourlyPmin][gT, 1])≈30
        @test normalized_rhs(EP1[:cHourlyBudget][1, gH])≈720
    end
end

end # module
