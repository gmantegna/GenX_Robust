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

    @testset "W1 two periods of unequal weight, hourly-varying Pmin and hydro maximum" begin
        # Review of GenX_Robust PR 14, finding 2: every other test has omega = 1 and profiles
        # constant within a period. Here T = 4, two periods of H = 2 h, Sub_Weights 4 and 12, so
        # omega = 4/2 = 2 in hours 1-2 and 12/2 = 6 in hours 3-4.
        # Zone 1: G (VRE 100 MW, cost 0, availability 1), demand 0. Zone 2: TH (thermal 100 MW,
        # 10 $/MWh, availability 1, 1, 1, 0.25), Hy (hydro 100 MW, budget column), demand 100 MW.
        # Line 1 -> 2: 100 MW, Line_Loss_Percentage 0.001 (half the loss at each end, so zone 2
        # receives 0.9995 of the sent flow), hurdle 5 $/MWh forward, 7 reverse.
        net_hdr = ",Network_zones,Network_Lines,Start_Zone,End_Zone,Line_Max_Flow_MW," *
                  "transmission_path_name,distance_mile,Line_Loss_Percentage," *
                  "Line_Max_Reinforcement_MW,Line_Reinforcement_Cost_per_MWyr," *
                  "Hurdle_Rate_Forward,Hurdle_Rate_Reverse"
        net_rows = [csvrow(["z1", "z1", 1, 1, 2, 100, "L12", 0, 0.001, 0, 0, 5.0, 7.0]),
            "z2,z2" * repeat(",", 11)]
        P = 2_000.0
        dir = ops_case(; T = 4, periods = 2, weights = [4, 12],
            demand = [[0.0, 0.0, 0.0, 0.0], [100.0, 100.0, 100.0, 100.0]],
            vre = [vre_row("G", 1, 100)], thermal = [thermal_row("TH", 2, 100, 10)],
            hydro = [hydro_row("Hy", 2, 100)],
            variability = Dict("TH" => [1.0, 1.0, 1.0, 0.25], "Hy" => [0.5, 0.2, 0.4, 0.1]),
            pmin = Dict("TH" => [0.1, 0.2, 0.3, 0.4]),
            budget = Dict("Hy" => [0.3, 0.3, 0.25, 0.25]),
            network = (net_hdr, net_rows))
        EP, inputs = ops_solve(dir; overrides = Dict("Hourly_Pmin" => 1,
            "HydroBudgetMaxProfile" => 1, "OperationalFloorPenalty" => P))
        @test termination_status(EP) == MOI.OPTIMAL
        @test inputs["omega"] == [2.0, 2.0, 6.0, 6.0]
        gT = findfirst(==("TH"), inputs["RESOURCE_NAMES"])
        gH = findfirst(==("Hy"), inputs["RESOURCE_NAMES"])
        obj = objective_function(EP)

        # Pmin rows: rhs = pi_t x 100 MW; hours 2, 3, 4: 0.2 x 100 = 20, 0.3 x 100 = 30,
        # 0.4 x 100 = 40 (hour 1: 0.1 x 100 = 10)
        @test normalized_rhs(EP[:cHourlyPmin][gT, 1])≈10 atol=ATOL
        @test normalized_rhs(EP[:cHourlyPmin][gT, 2])≈20 atol=ATOL
        @test normalized_rhs(EP[:cHourlyPmin][gT, 3])≈30 atol=ATOL
        @test normalized_rhs(EP[:cHourlyPmin][gT, 4])≈40 atol=ATOL
        # Hydro maximum rows: rhs = rho_t x 100 MW; hours 2, 3, 4: 0.2 x 100 = 20,
        # 0.4 x 100 = 40, 0.1 x 100 = 10 (hour 1: 0.5 x 100 = 50)
        @test normalized_rhs(EP[:cHydroBudgetMaxProfile][gH, 2])≈20 atol=ATOL
        @test normalized_rhs(EP[:cHydroBudgetMaxProfile][gH, 3])≈40 atol=ATOL
        @test normalized_rhs(EP[:cHydroBudgetMaxProfile][gH, 4])≈10 atol=ATOL

        # Objective coefficients. Floor slack: P x omega_t = 2,000 x 2 = 4,000 in hour 2,
        # 2,000 x 6 = 12,000 in hour 4. Hurdle: rate x omega_t = 5 x 2 = 10 forward in hour 2,
        # 5 x 6 = 30 forward and 7 x 6 = 42 reverse in hour 3.
        @test coefficient(obj, EP[:vFloorSlack][gT, 2])≈4_000 atol=ATOL
        @test coefficient(obj, EP[:vFloorSlack][gT, 4])≈12_000 atol=ATOL
        @test coefficient(obj, EP[:vTAUX_POS][1, 2])≈10 atol=ATOL
        @test coefficient(obj, EP[:vTAUX_POS][1, 3])≈30 atol=ATOL
        @test coefficient(obj, EP[:vTAUX_NEG][1, 3])≈42 atol=ATOL

        # Solution. Import costs 5 / 0.9995 = 5.0025 $/MWh delivered < 10 (TH), so TH runs at
        # its floor where it can: TH = 10, 20, 30 in hours 1-3; hour 4: floor 0.4 x 100 = 40,
        # availability 0.25 x 100 = 25, so TH = 25 and the floor slack u = 40 - 25 = 15.
        @test vp(EP, inputs, "TH")≈[10, 20, 30, 25] atol=ATOL
        @test value(EP[:vFloorSlack][gT, 4])≈15 atol=ATOL
        @test sum(value.(EP[:vFloorSlack]))≈15 atol=ATOL
        # Floor-slack cost: P x omega_4 x u = 2,000 x 6 x 15 = 180,000
        # (without the weight: 2,000 x 15 = 30,000)
        @test value(EP[:eTotalCFloorSlack])≈180_000 atol=1e-4
        # Hydro, period 2: budget (0.25 + 0.25) x 100 = 50 MWh = the maxima 40 + 10, so
        # Hy = 40 in hour 3 and 10 in hour 4. Period 1: budget (0.3 + 0.3) x 100 = 60 MWh.
        @test vp(EP, inputs, "Hy")[3:4]≈[40, 10] atol=ATOL
        @test sum(vp(EP, inputs, "Hy")[1:2])≈60 atol=ATOL
        @test sum(value.(EP[:vBudgetSpill]))≈0 atol=ATOL
        @test nse_mwh(EP)≈0 atol=ATOL
        # Delivered import = demand - TH - Hy: period 1: 200 - 30 - 60 = 110 MWh; hour 3:
        # 100 - 30 - 40 = 30; hour 4: 100 - 25 - 10 = 65. Weighted: 2 x 110 + 6 x (30 + 65) = 790.
        # Hurdle cost: 5 x 790 / 0.9995 = 3,951.97598799 (without the weight:
        # 5 x 205 / 0.9995 = 1,025.51275638)
        @test value(EP[:eTotalCHurdle])≈5 * 790 / 0.9995 atol=1e-6
        # Objective: TH 10 x (2 x (10 + 20) + 6 x (30 + 25)) = 3,900, + floor slack 180,000,
        # + hurdle 3,951.97598799 = 187,851.97598799
        @test objective_value(EP)≈3_900 + 180_000 + 5 * 790 / 0.9995 atol=1e-4
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
