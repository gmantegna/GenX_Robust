module TestOpsHurdle

# Hurdle rates on loss lines (Network.csv Hurdle_Rate_Forward/Hurdle_Rate_Reverse; reimagining-ra
# #204, spec mwe/specs/base_ops_port.md, PR item 4). Two zones, one line 1 -> 2 of 100 MW with
# Line_Loss_Percentage 0.001 (half the loss at each end), 1 hour, demand 50 MW.

include(joinpath(@__DIR__, "ops_rows_helpers.jl"))

const NET_HDR = ",Network_zones,Network_Lines,Start_Zone,End_Zone,Line_Max_Flow_MW,transmission_path_name,distance_mile,Line_Loss_Percentage,Line_Max_Reinforcement_MW,Line_Reinforcement_Cost_per_MWyr"

function net(; loss = 0.001, hf = 5.0, hr = 7.0, cols = (:f, :r))
    hdr = NET_HDR
    extra1 = String[]
    :f in cols && (hdr *= ",Hurdle_Rate_Forward"; push!(extra1, string(hf)))
    :r in cols && (hdr *= ",Hurdle_Rate_Reverse"; push!(extra1, string(hr)))
    r1 = csvrow(["z1", "z1", 1, 1, 2, 100, "L12", 0, loss, 0, 0])
    r2 = "z2,z2" * repeat(",", 9)
    if !isempty(extra1)
        r1 *= "," * join(extra1, ",")
        r2 *= repeat(",", length(extra1))
    end
    return (hdr, [r1, r2])
end

# G (100 MW, cost 0) in zone `gz`, demand 50 MW in the other zone
two_zone(gz; settings = Dict(), kw...) = ops_case(; T = 1,
    demand = gz == 1 ? [[0.0], [50.0]] : [[50.0], [0.0]],
    vre = [vre_row("G", gz, 100)], network = net(; kw...), settings = settings)

const PHI = 50 / (1 - 0.0005)   # 50.0250125 MW sent

@testset "Operational rows: hurdle rates" begin
    @testset "T6 forward: flow 50.0250125, hurdle 250.1250625" begin
        EP, _ = ops_solve(two_zone(1))
        @test termination_status(EP) == MOI.OPTIMAL
        @test value(EP[:vFLOW][1, 1])≈PHI atol=1e-7
        @test value(EP[:eTotalCHurdle])≈250.1250625 atol=1e-6
        @test objective_value(EP)≈250.1250625 atol=1e-6
    end

    @testset "T6 reverse: hurdle 350.1750875 on the end->start flow" begin
        EP, _ = ops_solve(two_zone(2))
        @test value(EP[:vFLOW][1, 1])≈-PHI atol=1e-7
        @test value(EP[:eTotalCHurdle])≈350.1750875 atol=1e-6
        @test objective_value(EP)≈350.1750875 atol=1e-6
    end

    @testset "T6 absent columns: no hurdle term, objective 0" begin
        EP, _ = ops_solve(two_zone(1; cols = ()))
        @test !has_name(EP, :eTotalCHurdle)
        @test objective_value(EP)≈0 atol=1e-6
    end

    @testset "T6 refusals, each with its accepted side" begin
        err(dir; o = Dict()) = errmsg_ops(() -> ops_inputs(dir, ops_setup(dir; overrides = o)))
        @test occursin("must be given together", err(two_zone(1; cols = (:f,))))
        @test err(two_zone(1; cols = (:f, :r))) == ""
        @test occursin("not a loss line", err(two_zone(1; loss = 0.0)))
        @test err(two_zone(1; loss = 0.0, hf = 0.0, hr = 0.0)) == ""
        @test occursin("Trans_Loss_Segments = 1 only",
            err(two_zone(1); o = Dict("Trans_Loss_Segments" => 0)))
        @test occursin("ParameterScale: 1", err(two_zone(1); o = Dict("ParameterScale" => 1)))
    end
end

end # module
