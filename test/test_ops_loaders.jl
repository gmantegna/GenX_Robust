module TestOpsLoaders

# Loaders, keys, refusals and Benders slicing of the operational rows ported from the chair's fork
# (reimagining-ra #204, spec mwe/specs/base_ops_port.md, PR item 1). Tiny cases written by
# ops_case (ops_rows_helpers.jl). No solve in this file.

include(joinpath(@__DIR__, "ops_rows_helpers.jl"))

const NAMES1 = ["G", "TH"]
base1(; kw...) = ops_case(; T = 1, demand = [50.0], vre = [vre_row("G", 1, 100)],
    thermal = [thermal_row("TH", 1, 100, 10)], kw...)

load_err(dir; overrides = Dict()) =
    errmsg_ops(() -> ops_inputs(dir, ops_setup(dir; overrides = overrides)))

@testset "Operational rows: loaders, keys, refusals, slicing" begin
    @testset "L1 absent inputs: nothing stored, defaults" begin
        dir = base1()
        setup = ops_setup(dir)
        inputs = ops_inputs(dir, setup)
        @test !haskey(inputs, "pP_Min")
        @test !haskey(inputs, "pBudget")
        d = GenX.default_settings()
        @test d["Hourly_Pmin"] == 0
        @test d["Fixed_Dispatch"] == 0
        @test d["HydroBudgetMaxProfile"] == 0
        @test d["StorageDischargeDerate"] == 0
        @test d["StorageChargeHeadroom"] == 1
        @test !haskey(d, "OperationalFloorPenalty")
        @test !haskey(d, "EnergyBudgetWindowHours")
    end

    @testset "L2 Generators_Pmin.csv read with Hourly_Pmin 1; missing column = 0" begin
        dir = base1(pmin = Dict("TH" => [0.3]))
        inputs = ops_inputs(dir, ops_setup(dir; overrides = Dict("Hourly_Pmin" => 1)))
        gT = findfirst(==("TH"), inputs["RESOURCE_NAMES"])
        gG = findfirst(==("G"), inputs["RESOURCE_NAMES"])
        @test inputs["pP_Min"][gT, 1] == 0.3
        @test inputs["pP_Min"][gG, 1] == 0.0
        @test inputs["PMIN_RIDS"] == [gT]
        # Hourly_Pmin 0: the file is not read
        inputs0 = ops_inputs(dir, ops_setup(dir))
        @test !haskey(inputs0, "pP_Min")
    end

    @testset "L3 Hourly_Pmin 1 needs the file" begin
        dir = base1()
        @test occursin("Generators_Pmin.csv does not exist",
            load_err(dir; overrides = Dict("Hourly_Pmin" => 1)))
        @test load_err(dir) == ""
    end

    @testset "L4 Hourly_energy_budget.csv read when present" begin
        dir = base1(budget = Dict("G" => [0.5]))
        inputs = ops_inputs(dir, ops_setup(dir))
        gG = findfirst(==("G"), inputs["RESOURCE_NAMES"])
        gT = findfirst(==("TH"), inputs["RESOURCE_NAMES"])
        @test inputs["pBudget"][gG, 1] == 0.5
        @test inputs["pBudget"][gT, 1] == 0.0
        @test inputs["BUDGET_RIDS"] == [gG]
    end

    @testset "L5 Fixed_Dispatch 1 refused, 0 accepted" begin
        dir = base1()
        @test occursin("Fixed_Dispatch: 1",
            load_err(dir; overrides = Dict("Fixed_Dispatch" => 1)))
        @test load_err(dir; overrides = Dict("Fixed_Dispatch" => 0)) == ""
        @test occursin("must be 0 or 1",
            load_err(dir; overrides = Dict("HydroBudgetMaxProfile" => 2)))
    end

    @testset "L6 ParameterScale 1 refused with these inputs, accepted without" begin
        dir = base1(pmin = Dict("TH" => [0.3]))
        @test occursin("ParameterScale: 1 is not supported",
            load_err(dir; overrides = Dict("Hourly_Pmin" => 1, "ParameterScale" => 1)))
        @test load_err(dir; overrides = Dict("Hourly_Pmin" => 1, "ParameterScale" => 0)) == ""
        @test load_err(dir; overrides = Dict("Hourly_Pmin" => 0, "ParameterScale" => 1)) == ""
        dirb = base1(budget = Dict("G" => [0.5]))
        @test occursin("ParameterScale: 1 is not supported",
            load_err(dirb; overrides = Dict("ParameterScale" => 1)))
    end

    @testset "L7 OperationalFloorPenalty must exceed VoLL (1000)" begin
        dir = base1()
        @test occursin("must be a number above VoLL",
            load_err(dir; overrides = Dict("OperationalFloorPenalty" => 1000)))
        @test load_err(dir; overrides = Dict("OperationalFloorPenalty" => 1000.5)) == ""
        @test load_err(dir; overrides = Dict("OperationalFloorPenalty" => 1001)) == ""
        @test occursin("must be a number above VoLL",
            load_err(dir; overrides = Dict("OperationalFloorPenalty" => "high")))
    end

    @testset "L7b VoLL is the largest segment cost (two segments: 1000 and 1500)" begin
        dir = ops_case(; T = 2, demand = [50.0, 50.0], vre = [vre_row("G", 1, 100)])
        f = joinpath(dir, "system", "Demand_data.csv")
        lines = readlines(f)
        c = split(lines[3], ",")
        c[2] = "2"; c[3] = "1.5"; c[4] = "0.1"
        lines[3] = join(c, ",")
        write(f, join(lines, "\n") * "\n")
        @test occursin("above VoLL = 1500",
            load_err(dir; overrides = Dict("OperationalFloorPenalty" => 1200)))
        @test load_err(dir; overrides = Dict("OperationalFloorPenalty" => 1600)) == ""
    end

    @testset "L8 EnergyBudgetWindowHours divides the subperiod (H = 48)" begin
        dir = ops_case(; T = 48, demand = fill(50.0, 48), vre = [vre_row("G", 1, 100)])
        @test load_err(dir; overrides = Dict("EnergyBudgetWindowHours" => 24)) == ""
        @test load_err(dir; overrides = Dict("EnergyBudgetWindowHours" => 48)) == ""
        @test occursin("EnergyBudgetWindowHours must be",
            load_err(dir; overrides = Dict("EnergyBudgetWindowHours" => 25)))
        @test occursin("EnergyBudgetWindowHours must be",
            load_err(dir; overrides = Dict("EnergyBudgetWindowHours" => 0)))
        @test occursin("EnergyBudgetWindowHours must be",
            load_err(dir; overrides = Dict("EnergyBudgetWindowHours" => 96)))
    end

    @testset "L9 budgeted hydro with LDS 1 refused, LDS 0 accepted" begin
        mk(lds) = ops_case(; T = 24, demand = fill(60.0, 24),
            hydro = [hydro_row("Hy", 1, 100; lds = lds)], budget = Dict("Hy" => fill(0.5, 24)))
        @test occursin("cannot have LDS = 1", load_err(mk(1)))
        @test load_err(mk(0)) == ""
    end

    @testset "L10 OperationalFloorPenalty refused with NSEBudget 1, accepted with NSEBudget 0" begin
        # The penalty's guard (P > VoLL) assumes shed is priced at VoLL; under the NSE budget it
        # is not (review of GenX_Robust PR 14, finding 1).
        dir = base1()
        nse1 = Dict("NSEBudget" => 1, "NSEBudgetTargetMWh" => 0.0)
        @test occursin("NSE-budget runs must use hard floors and budgets",
            load_err(dir; overrides = merge(nse1, Dict("OperationalFloorPenalty" => 2_000))))
        @test load_err(dir;
            overrides = Dict("NSEBudget" => 0, "OperationalFloorPenalty" => 2_000)) == ""
        @test load_err(dir; overrides = nse1) == ""
    end

    @testset "T11 slicing of pP_Min and pBudget into Benders subperiods" begin
        # 2 periods of 24 h, capacity 100 MW: Pmin 0.3 / 0.6, budget 0.3 / 0.6 per hour
        prof = [fill(0.3, 24); fill(0.6, 24)]
        dir = ops_case(; T = 48, periods = 2, demand = fill(50.0, 48),
            thermal = [thermal_row("TH", 1, 100, 10)], hydro = [hydro_row("Hy", 1, 100)],
            pmin = Dict("TH" => prof), budget = Dict("Hy" => prof))
        inputs = ops_inputs(dir, ops_setup(dir; overrides = Dict("Hourly_Pmin" => 1)))
        gT = findfirst(==("TH"), inputs["RESOURCE_NAMES"])
        gH = findfirst(==("Hy"), inputs["RESOURCE_NAMES"])
        parts = GenX.separate_inputs_subperiods(inputs)
        @test size(parts[2]["pP_Min"]) == (inputs["G"], 24)
        @test all(parts[2]["pP_Min"][gT, :] .* 100 .≈ 60)   # unsliced: 30
        @test all(parts[1]["pP_Min"][gT, :] .* 100 .≈ 30)
        @test sum(parts[2]["pBudget"][gH, :]) * 100 ≈ 1440  # unsliced: 720
        @test sum(parts[1]["pBudget"][gH, :]) * 100 ≈ 720
    end
end

end # module
