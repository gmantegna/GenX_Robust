module TestOpsSimCharge

# Simultaneous charging (StorSimultaneousCharge, default 0; reimagining-ra #204): the row
# cStorSimCharge of the lab's rung-B branch base/d58-rungb f63283ac
# (src/model/resources/storage/storage_symmetric.jl there):
#   vP[y,t] + vCHARGE[y,t] <= 0.5 x (pP_Max[y,t] + 1) x capacity
# for operational symmetric storage whose pP_Max is below 1 in some hour (a group's representative
# reads the group capacity eTotalCapStorOp). Hand values in the comments; omega = 1, VoLL 1000.

include(joinpath(@__DIR__, "ops_rows_helpers.jl"))

const ATOL = 1e-6

# 3 hours, efficiency 1. G (100 MW, cost 0) produces only in hour 2; demand only in hour 3
# (`d3` MW). Storage pP_Max (Generators_variability.csv) is given per storage name in `svar`.
# StorageDischargeDerate stays 0 (default), so the profile acts only through cStorSimCharge.
sim_case(; storage, svar, d3, group = nothing) = ops_case(; T = 3, demand = [0.0, 0.0, d3],
    vre = [vre_row("G", 1, 100)], storage = storage, stor_group = group,
    variability = merge(Dict{String, Any}("G" => [0.0, 1.0, 0.0]), svar))

rid(inputs, n) = findfirst(==(n), inputs["RESOURCE_NAMES"])

const PROF = [1.0, 0.5, 1.0]   # below 1 in hour 2 only, 1 in hours 1 and 3

@testset "Operational rows: simultaneous charging (StorSimultaneousCharge)" begin
    @testset "S1 ungrouped A (10 MW, 40 MWh): the row binds in hour 2, not the first" begin
        dir = sim_case(storage = [stor_row("A", 1, 10, 40)], svar = Dict("A" => PROF), d3 = 20.0)
        # Key 0: charge in hour 2 <= 10 (vP + vCHARGE <= 10), discharge 10 in hour 3,
        # NSE 20 - 10 = 10 MWh, objective 10 x 1000 = 10,000.
        EP0, inputs = ops_solve(dir)
        @test termination_status(EP0) == MOI.OPTIMAL
        @test !has_name(EP0, :cStorSimCharge)
        @test objective_value(EP0)≈10_000 atol=ATOL
        # Key 1: hour 2 row vP + vCHARGE <= 0.5 x (0.5 + 1) x 10 = 7.5, so charge 7.5,
        # discharge 7.5 in hour 3, NSE 20 - 7.5 = 12.5, objective 12,500.
        EP1, _ = ops_solve(dir; overrides = Dict("StorSimultaneousCharge" => 1))
        @test termination_status(EP1) == MOI.OPTIMAL
        gA = rid(inputs, "A")
        @test normalized_rhs(EP1[:cStorSimCharge][gA, 2])≈7.5 atol=ATOL  # 0.5 x 1.5 x 10
        @test normalized_rhs(EP1[:cStorSimCharge][gA, 1])≈10 atol=ATOL   # 0.5 x 2 x 10
        @test normalized_rhs(EP1[:cStorSimCharge][gA, 3])≈10 atol=ATOL
        @test value(EP1[:vCHARGE][gA, 2])≈7.5 atol=ATOL
        @test nse_mwh(EP1)≈12.5 atol=ATOL
        @test objective_value(EP1)≈12_500 atol=ATOL
    end

    @testset "S2 profile 1 in every hour: no row, objective as with the key off" begin
        dir = sim_case(storage = [stor_row("A", 1, 10, 40)], svar = Dict("A" => ones(3)), d3 = 20.0)
        EP1, _ = ops_solve(dir; overrides = Dict("StorSimultaneousCharge" => 1))
        @test isempty(EP1[:cStorSimCharge])
        @test objective_value(EP1)≈10_000 atol=ATOL      # as S1 key 0
    end

    @testset "S3 group S = A (10 MW, 40 MWh, representative) + B (30 MW, 120 MWh)" begin
        st = [stor_row("A", 1, 10, 40), stor_row("B", 1, 30, 120)]
        # B's column (0.5 every hour) is a member's column: not read.
        dir = sim_case(storage = st, group = ["S", "S"], d3 = 40.0,
            svar = Dict("A" => PROF, "B" => fill(0.5, 3)))
        # Key 0: group capacity 10 + 30 = 40: charge 40 in hour 2, discharge 40, NSE 0, objective 0.
        EP0, inputs = ops_solve(dir)
        @test objective_value(EP0)≈0 atol=ATOL
        # Key 1: hour 2 row on A: vP + vCHARGE <= 0.5 x 1.5 x 40 = 30; NSE 40 - 30 = 10,
        # objective 10,000. (On A's own 10 MW: 7.5, NSE 32.5, objective 32,500.)
        EP1, _ = ops_solve(dir; overrides = Dict("StorSimultaneousCharge" => 1))
        gA, gB = rid(inputs, "A"), rid(inputs, "B")
        @test collect(axes(EP1[:cStorSimCharge], 1)) == [gA]   # the member B has no row
        @test normalized_rhs(EP1[:cStorSimCharge][gA, 2])≈30 atol=ATOL
        @test nse_mwh(EP1)≈10 atol=ATOL
        @test objective_value(EP1)≈10_000 atol=ATOL
    end

    @testset "S4 the same A and B in no group: each carries its own row" begin
        st = [stor_row("A", 1, 10, 40), stor_row("B", 1, 30, 120)]
        dir = sim_case(storage = st, group = ["None", "None"], d3 = 40.0,
            svar = Dict("A" => PROF, "B" => ones(3)))
        # Key 1: A charges <= 0.5 x 1.5 x 10 = 7.5 in hour 2; B (profile 1, no row) charges 30;
        # discharge 37.5, NSE 2.5, objective 2,500.
        EP1, inputs = ops_solve(dir; overrides = Dict("StorSimultaneousCharge" => 1))
        @test collect(axes(EP1[:cStorSimCharge], 1)) == [rid(inputs, "A")]
        @test objective_value(EP1)≈2_500 atol=ATOL
    end

    @testset "S5 Benders subproblem builds cStorSimCharge (same values as the monolith)" begin
        dir = sim_case(storage = [stor_row("A", 1, 10, 40)], svar = Dict("A" => PROF), d3 = 20.0)
        setup = ops_setup(dir; overrides = Dict("StorSimultaneousCharge" => 1, "Benders" => 1))
        inputs = ops_inputs(dir, setup)
        parts = GenX.separate_inputs_subperiods(inputs)
        OPT = optimizer_with_attributes(HiGHS.Optimizer)
        EP = quiet_ops(() -> GenX.generate_operation_subproblem(setup, parts[1], OPT))
        set_silent(EP)
        gA = rid(inputs, "A")
        @test normalized_rhs(EP[:cStorSimCharge][gA, 2])≈7.5 atol=ATOL
        optimize!(EP)
        @test objective_value(EP)≈12_500 atol=ATOL       # as S1 key 1
    end

    @testset "S6 the key must be 0 or 1" begin
        dir = sim_case(storage = [stor_row("A", 1, 10, 40)], svar = Dict("A" => PROF), d3 = 20.0)
        err(o) = errmsg_ops(() -> ops_inputs(dir, ops_setup(dir; overrides = o)))
        @test occursin("StorSimultaneousCharge must be 0 or 1",
            err(Dict("StorSimultaneousCharge" => 2)))
        @test err(Dict("StorSimultaneousCharge" => 1)) == ""
    end
end

end # module
