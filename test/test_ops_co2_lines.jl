module TestOpsCO2Lines

# Import CO2 on lines for the mass cap (Network.csv EF_Forward_CO_2_Cap_Zone_c /
# EF_Reverse_CO_2_Cap_Zone_c; ported from the chair's fork d58adffd co2_cap.jl:93-130;
# reimagining-ra #204, spec mwe/specs/base_ops_port.md F7f and T7). Two zones, one line of 100 MW
# with Line_Loss_Percentage 0.001 (half the loss at each end), 1 hour, omega = 1. G (100 MW, cost 0,
# no emissions) in zone 1; TH (100 MW, 10 $/MWh, no emissions) in zone 2; demand 50 MW in zone 2.
# One CO2 cap over both zones of 10 t (CO_2_Max_Mtons 1e-5). Emissions come only from the import
# term, so the cap binds only because of it.

include(joinpath(@__DIR__, "ops_rows_helpers.jl"))

const ATOL = 1e-6
const NET_HDR = ",Network_zones,Network_Lines,Start_Zone,End_Zone,Line_Max_Flow_MW,transmission_path_name,distance_mile,Line_Loss_Percentage,Line_Max_Reinforcement_MW,Line_Reinforcement_Cost_per_MWyr"

# line start -> stop; factors ef (forward), er (reverse); cols = which prefixes are written
function net(; start = 1, stop = 2, loss = 0.001, ef = 0.5, er = 0.0, cols = (:f, :r))
    hdr = NET_HDR
    extra = String[]
    :f in cols && (hdr *= ",EF_Forward_CO_2_Cap_Zone_1"; push!(extra, string(ef)))
    :r in cols && (hdr *= ",EF_Reverse_CO_2_Cap_Zone_1"; push!(extra, string(er)))
    r1 = csvrow(["z1", "z1", 1, start, stop, 100, "L", 0, loss, 0, 0])
    r2 = "z2,z2" * repeat(",", 9)
    if !isempty(extra)
        r1 *= "," * join(extra, ",")
        r2 *= repeat(",", length(extra))
    end
    return (hdr, [r1, r2])
end

function co2_case(; settings = Dict("CO2Cap" => 1), weights = nothing, kw...)
    dir = ops_case(; T = 1, demand = [[0.0], [50.0]], vre = [vre_row("G", 1, 100)], weights = weights,
        thermal = [thermal_row("TH", 2, 100, 10)], network = net(; kw...), settings = settings)
    write(joinpath(dir, "policies", "CO2_cap.csv"),
        ",Network_zones,CO_2_Cap_Zone_1,CO_2_Max_tons_MWh_1,CO_2_Max_Mtons_1\n" *
        "z1,z1,1,0.05,1e-5\nz2,z2,1,0.05,0\n")
    return dir
end

err(dir; o = Dict()) = errmsg_ops(() -> ops_inputs(dir, ops_setup(dir; overrides = o)))

@testset "Operational rows: import CO2 on lines" begin
    @testset "T7 forward: Phi+ = 20, TH = 30.01, objective 300.1" begin
        EP, _ = ops_solve(co2_case())
        @test termination_status(EP) == MOI.OPTIMAL
        # 0.5 t/MWh x Phi+ <= 10 t -> Phi+ = 20 MW sent; zone 2 receives 20 - 0.0005 x 20 = 19.99;
        # TH = 50 - 19.99 = 30.01 MW; objective 10 x 30.01 = 300.1
        @test value(EP[:vTAUX_POS][1, 1])≈20 atol=ATOL
        @test value(EP[:eCO2LineImports][1])≈10 atol=ATOL
        @test value(EP[:eForwardEmissionsByLine][1, 1])≈10 atol=ATOL
        @test value(EP[:cCO2Emissions_systemwide][1])≈10 atol=ATOL       # the row binds
        @test objective_value(EP)≈300.1 atol=ATOL
        # coefficient of vTAUX_POS[1,1] in the row = omega x EF+ = 1 x 0.5
        @test normalized_coefficient(EP[:cCO2Emissions_systemwide][1], EP[:vTAUX_POS][1, 1])≈0.5
    end

    @testset "T7 weighted hour (omega = 2): Phi+ = 10, objective 800.1" begin
        # one 1-hour period with Sub_Weights 2: omega = 2 / 1 = 2
        dir = co2_case(; weights = [2])
        EP, inputs = ops_solve(dir)
        @test inputs["omega"] == [2.0]
        # 2 x 0.5 x Phi+ <= 10 t -> Phi+ = 10; zone 2 receives 10 x 0.9995 = 9.995; TH = 40.005;
        # objective 2 x 10 x 40.005 = 800.1 (without omega in the term: Phi+ = 20, 600.2)
        @test value(EP[:vTAUX_POS][1, 1])≈10 atol=ATOL
        @test objective_value(EP)≈800.1 atol=ATOL
        @test normalized_coefficient(EP[:cCO2Emissions_systemwide][1], EP[:vTAUX_POS][1, 1])≈1.0
    end

    @testset "T7 reverse: line 2 -> 1, flow is Phi-, EF+ 0.1, EF- 0.5: the same values" begin
        EP, _ = ops_solve(co2_case(; start = 2, stop = 1, ef = 0.1, er = 0.5))
        # 0.1 x Phi+ + 0.5 x Phi- <= 10 with Phi+ = 0: Phi- = 20, TH = 30.01, objective 300.1
        # (EF+ applied to Phi- would allow 100 MW and give objective 0)
        @test value(EP[:vTAUX_NEG][1, 1])≈20 atol=ATOL
        @test value(EP[:eReverseEmissionsByLine][1, 1])≈10 atol=ATOL
        @test objective_value(EP)≈300.1 atol=ATOL
    end

    @testset "T7 S4: a reverse factor with a zero forward factor adds nothing (fork's filter)" begin
        EP, _ = ops_solve(co2_case(; start = 2, stop = 1, ef = 0.0, er = 0.5))
        # line filtered out of the cap: import 50 / 0.9995 = 50.025 MW, no TH, objective 0
        @test value(EP[:eCO2LineImports][1])≈0 atol=ATOL
        @test objective_value(EP)≈0 atol=ATOL
    end

    @testset "T7 absent columns and all-zero factors: no term, objective 0" begin
        EP, _ = ops_solve(co2_case(; cols = ()))
        @test !has_name(EP, :eCO2LineImports)
        @test objective_value(EP)≈0 atol=ATOL
        EP0, _ = ops_solve(co2_case(; ef = 0.0, er = 0.0))
        @test !has_name(EP0, :eCO2LineImports)
        @test objective_value(EP0)≈0 atol=ATOL
    end

    @testset "T7 Benders subproblem: the row carries the term, budget 10 t gives 300.1" begin
        dir = co2_case()
        setup = ops_setup(dir; overrides = Dict("Benders" => 1))
        inputs = ops_inputs(dir, setup)
        parts = GenX.separate_inputs_subperiods(inputs)
        OPT = optimizer_with_attributes(HiGHS.Optimizer)
        EP = quiet_ops(() -> GenX.generate_operation_subproblem(setup, parts[1], OPT))
        set_silent(EP)
        c = EP[:cCO2Emissions_systemwide][1]
        @test normalized_coefficient(c, EP[:vTAUX_POS][1, 1])≈0.5
        fix(EP[:vCO2budget][1, 1], 10.0; force = true)       # the master's budget for period 1
        optimize!(EP)
        @test objective_value(EP)≈300.1 atol=ATOL
    end

    @testset "T7 refusals, each with its accepted side" begin
        @test occursin("must be given together", err(co2_case(; cols = (:f,))))
        @test err(co2_case()) == ""
        @test occursin("not a loss line", err(co2_case(; loss = 0.0)))
        @test err(co2_case(; loss = 0.0, ef = 0.0, er = 0.0)) == ""
        @test occursin("Trans_Loss_Segments = 1 only",
            err(co2_case(); o = Dict("Trans_Loss_Segments" => 0)))
        @test occursin("not supported with ParameterScale: 1",
            err(co2_case(); o = Dict("ParameterScale" => 1)))
        @test occursin("mass cap (CO2Cap: 1) only", err(co2_case(); o = Dict("CO2Cap" => 2)))
        @test err(co2_case(; ef = 0.0, er = 0.0); o = Dict("CO2Cap" => 2)) == ""
    end
end

end # module
