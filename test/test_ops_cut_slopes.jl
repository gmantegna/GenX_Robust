module TestOpsCutSlopes

# Benders cut slopes of the ported operational rows (reimagining-ra #204; review of GenX PR 14,
# finding "Benders cut coefficients of the new rows are untested"). For a buildable resource the
# subproblem's linking variable (vCAP[g], vCO2budget[w,c]) is fixed as MacroEnergySolvers does
# (init_subproblem, then fix(...; force = true)); the cut coefficient is dual(FixRef(v)), the value
# MES puts in the cut (src/benders/mes_feasible_point_shim.jl). It must equal the finite
# difference of the subproblem's optimal value V. The step h is checked not to cross a basis
# change: the forward difference (V(x+h)-V(x))/h and the backward difference (V(x)-V(x-h))/h are
# both asserted equal to the dual. V is convex and piecewise linear in x, so equal one-sided
# differences mean V is linear on [x-h, x+h]. Hand values in the comments; omega = 1, VoLL 1000.

include(joinpath(@__DIR__, "ops_rows_helpers.jl"))

const ATOL = 1e-6

"Set column `col` of the CSV row `row` written for header `hdr` to `val`."
function setcol(row, hdr, col, val)
    h = split(hdr, ",")
    r = split(row, ",")
    r[findfirst(==(col), h)] = string(val)
    return join(r, ",")
end

rid(inputs, n) = findfirst(==(n), inputs["RESOURCE_NAMES"])

"""
Subproblem of period 1 with the linking variables `fixes` (name => value) fixed; returns the
model after optimize!.
"""
function solve_sub(dir, overrides, fixes)
    setup = ops_setup(dir; overrides = merge(Dict{String, Any}("Benders" => 1), overrides))
    inputs = ops_inputs(dir, setup)
    parts = GenX.separate_inputs_subperiods(inputs)
    OPT = optimizer_with_attributes(HiGHS.Optimizer)
    EP, linking = quiet_ops(() -> GenX.init_subproblem(setup, parts[1], OPT, collect(keys(fixes))))
    sort(linking) == sort(collect(keys(fixes))) || error("linking variables $linking")
    for (n, v) in fixes
        fix(variable_by_name(EP, n), v; force = true)
    end
    set_attribute(EP, "threads", 1)
    optimize!(EP)
    return EP
end

"""
(dual, forward difference, backward difference) of V in the linking variable `var` at `x`, the
other linking variables at `others`.
"""
function slope(dir, overrides, var, x, h; others = Dict{String, Float64}())
    EP = solve_sub(dir, overrides, merge(others, Dict(var => x)))
    termination_status(EP) == MOI.OPTIMAL || error("subproblem not optimal")
    V0 = objective_value(EP)
    lam = dual(FixRef(variable_by_name(EP, var)))
    Vp = objective_value(solve_sub(dir, overrides, merge(others, Dict(var => x + h))))
    Vm = objective_value(solve_sub(dir, overrides, merge(others, Dict(var => x - h))))
    return (V0 = V0, lam = lam, fwd = (Vp - V0) / h, bwd = (V0 - Vm) / h)
end

buildable_thermal(n, z, varom) = setcol(thermal_row(n, z, 0, varom), THERMAL_HDR, "New_Build", 1)
buildable_hydro(n, z) = setcol(hydro_row(n, z, 0), HYDRO_HDR, "New_Build", 1)
buildable_stor(n, z) = setcol(stor_row(n, z, 0, 0; mind = 4, maxd = 4), STOR_HDR, "New_Build", 1)

@testset "Operational rows: Benders cut slopes of the new rows" begin
    @testset "C1 hourly floor on a buildable thermal" begin
        # 2 h, demand 50, G 100 MW at cost 0, TH buildable at 10 $/MWh with Pmin 0.3, 0.5.
        # At capacity D = 40 the floors (12, 20 MW) force TH: V = 10 x (0.3 + 0.5) x D = 320,
        # dV/dD = 10 x 0.8 = 8.
        dir = ops_case(; T = 2, demand = [50.0, 50.0], vre = [vre_row("G", 1, 100)],
            thermal = [buildable_thermal("TH", 1, 10)], pmin = Dict("TH" => [0.3, 0.5]))
        inputs = ops_inputs(dir, ops_setup(dir))
        s = slope(dir, Dict("Hourly_Pmin" => 1), "vCAP[$(rid(inputs, "TH"))]", 40.0, 0.5)
        @test s.V0≈320 atol=ATOL
        @test s.lam≈8 atol=ATOL
        @test s.fwd≈8 atol=ATOL
        @test s.bwd≈8 atol=ATOL
    end

    @testset "C2 energy budget on a buildable hydro" begin
        # 4 h, demand 50, TH 100 MW existing at 10 $/MWh, Hy buildable with budget column
        # 0.2, 0.4, 0.6, 0.8 (sum 2.0). At D = 20 Hy must produce 2.0 x 20 = 40 MWh (<= 20 MW
        # each hour): V = 10 x (200 - 2.0 D) = 1,600, dV/dD = -10 x 2.0 = -20.
        dir = ops_case(; T = 4, demand = fill(50.0, 4), thermal = [thermal_row("TH", 1, 100, 10)],
            hydro = [buildable_hydro("Hy", 1)], budget = Dict("Hy" => [0.2, 0.4, 0.6, 0.8]))
        inputs = ops_inputs(dir, ops_setup(dir))
        s = slope(dir, Dict(), "vCAP[$(rid(inputs, "Hy"))]", 20.0, 0.5)
        @test s.V0≈1_600 atol=ATOL
        @test s.lam≈-20 atol=ATOL
        @test s.fwd≈-20 atol=ATOL
        @test s.bwd≈-20 atol=ATOL
    end

    st_case(; T, demand, g, avar) = ops_case(; T = T, demand = demand,
        vre = [vre_row("G", 1, 100)],
        storage = [stor_row("A", 1, 10, 40), buildable_stor("B", 1)], stor_group = ["S", "S"],
        variability = Dict("G" => g, "A" => avar))

    @testset "C3 storage group: buildable member B, derate on the representative" begin
        # 2 h; G 100 MW in h1, 0 in h2; demand 0, 35. Group S = A (10 MW, 40 MWh, representative)
        # + B (buildable, power D, energy 80 MWh). StorageDischargeDerate 1, A's column 1, 0.5.
        # Discharge in h2 <= 0.5 x (10 + D); at D = 20: 15 MW, NSE 20, V = 20,000;
        # dV/dD = -1000 x 0.5 = -500.
        dir = st_case(T = 2, demand = [0.0, 35.0], g = [1.0, 0.0], avar = [1.0, 0.5])
        inputs = ops_inputs(dir, ops_setup(dir))
        gB = rid(inputs, "B")
        s = slope(dir, Dict("StorageDischargeDerate" => 1), "vCAP[$gB]", 20.0, 0.5;
            others = Dict("vCAPENERGY[$gB]" => 80.0))
        @test s.V0≈20_000 atol=ATOL
        @test s.lam≈-500 atol=ATOL
        @test s.fwd≈-500 atol=ATOL
        @test s.bwd≈-500 atol=ATOL
    end

    @testset "C4 simultaneous charging (B1) on the group: buildable member B" begin
        # 3 h; G 100 MW only in h2; demand 40 in h3 only. Group S = A (10 MW, 40 MWh) + B
        # (power D, energy 120 MWh); A's column 1, 0.5, 1; StorSimultaneousCharge 1, no derate.
        # Charge in h2 <= 0.5 x (0.5 + 1) x (10 + D) = 0.75 (10 + D); at D = 30: 30, NSE 10,
        # V = 10,000; dV/dD = -1000 x 0.75 = -750.
        dir = st_case(T = 3, demand = [0.0, 0.0, 40.0], g = [0.0, 1.0, 0.0], avar = [1.0, 0.5, 1.0])
        inputs = ops_inputs(dir, ops_setup(dir))
        gB = rid(inputs, "B")
        s = slope(dir, Dict("StorSimultaneousCharge" => 1), "vCAP[$gB]", 30.0, 0.5;
            others = Dict("vCAPENERGY[$gB]" => 120.0))
        @test s.V0≈10_000 atol=ATOL
        @test s.lam≈-750 atol=ATOL
        @test s.fwd≈-750 atol=ATOL
        @test s.bwd≈-750 atol=ATOL
    end

    @testset "C5 import CO2 on lines (B2): slope in the period's CO2 budget" begin
        # Two zones, line 1 -> 2 (100 MW, loss 0.001), EF+ 0.5 t/MWh; G (cost 0) in zone 1, TH
        # (10 $/MWh) and demand 50 in zone 2; CO2Cap 1. Budget b t: Phi+ = 2b MW sent,
        # 2b x 0.9995 delivered, TH = 50 - 1.999 b; at b = 10: V = 300.1;
        # dV/db = -10 x 2 x 0.9995 = -19.99.
        hdr = ",Network_zones,Network_Lines,Start_Zone,End_Zone,Line_Max_Flow_MW,transmission_path_name,distance_mile,Line_Loss_Percentage,Line_Max_Reinforcement_MW,Line_Reinforcement_Cost_per_MWyr,EF_Forward_CO_2_Cap_Zone_1,EF_Reverse_CO_2_Cap_Zone_1"
        rows = [csvrow(["z1", "z1", 1, 1, 2, 100, "L", 0, 0.001, 0, 0, 0.5, 0.0]),
            "z2,z2" * repeat(",", 11)]
        dir = ops_case(; T = 1, demand = [[0.0], [50.0]], vre = [vre_row("G", 1, 100)],
            thermal = [thermal_row("TH", 2, 100, 10)], network = (hdr, rows),
            settings = Dict("CO2Cap" => 1))
        write(joinpath(dir, "policies", "CO2_cap.csv"),
            ",Network_zones,CO_2_Cap_Zone_1,CO_2_Max_tons_MWh_1,CO_2_Max_Mtons_1\n" *
            "z1,z1,1,0.05,1e-5\nz2,z2,1,0.05,0\n")
        s = slope(dir, Dict(), "vCO2budget[1,1]", 10.0, 0.5)
        @test s.V0≈300.1 atol=ATOL
        @test s.lam≈-19.99 atol=ATOL
        @test s.fwd≈-19.99 atol=ATOL
        @test s.bwd≈-19.99 atol=ATOL
    end
end

end # module
