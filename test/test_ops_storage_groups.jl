module TestOpsStorageGroups

# Storage operational groups (Stor_Op_Group), the storage discharge derate (StorageDischargeDerate)
# and the charge-headroom switch (StorageChargeHeadroom) (reimagining-ra #204, spec
# mwe/specs/base_ops_port.md, PR item 3). Hand values from the spec's section 9.

include(joinpath(@__DIR__, "ops_rows_helpers.jl"))

const ATOL = 1e-6

# A (10 MW, 40 MWh) and B (30 MW, 120 MWh), efficiency 1, generator G 100 MW (profile g)
function ab_case(; demand, g, group = ["S", "S"], a_var = nothing, a_kw = (;), b_kw = (;),
        settings = Dict())
    var = Dict{String, Any}("G" => g)
    a_var === nothing || (var["A"] = a_var)
    return ops_case(; T = length(demand), demand = demand, vre = [vre_row("G", 1, 100)],
        storage = [stor_row("A", 1, 10, 40; a_kw...), stor_row("B", 1, 30, 120; b_kw...)],
        stor_group = group, variability = var, settings = settings)
end

ncons(EP) = num_constraints(EP; count_variable_in_set_constraints = false)

@testset "Operational rows: storage groups, derate, charge headroom" begin
    @testset "T5 group + derate: 0 without derate, NSE 15 / 15,000 with it" begin
        dir = ab_case(demand = [0.0, 35.0], g = [1.0, 0.0], a_var = [1.0, 0.5])
        EP0, inputs = ops_solve(dir; overrides = Dict("StorageDischargeDerate" => 0))
        @test termination_status(EP0) == MOI.OPTIMAL
        @test objective_value(EP0)≈0 atol=ATOL
        @test nse_mwh(EP0)≈0 atol=ATOL
        EP1, _ = ops_solve(dir; overrides = Dict("StorageDischargeDerate" => 1))
        @test nse_mwh(EP1)≈15 atol=ATOL
        @test objective_value(EP1)≈15_000 atol=ATOL
        @test vp(EP1, inputs, "A")[2]≈20 atol=ATOL
        @test all(abs.(vp(EP1, inputs, "B")) .<= ATOL)          # Theta_B = 0
        gB = findfirst(==("B"), inputs["RESOURCE_NAMES"])
        @test all(abs.(value.(EP1[:vCHARGE][gB, :])) .<= ATOL)
        @test inputs["STOR_BIN_ONLY"] == [gB]
        @test value(EP1[:eTotalCapStorOp][findfirst(==("A"), inputs["RESOURCE_NAMES"])])≈40
        @test value(EP1[:eTotalCapEnergyStorOp][findfirst(==("A"), inputs["RESOURCE_NAMES"])])≈160
    end

    @testset "T5 ungrouped storage gets the derate too" begin
        dir = ab_case(demand = [0.0, 35.0], g = [1.0, 0.0], a_var = [1.0, 0.5],
            group = ["None", "None"])
        EP1, inputs = ops_solve(dir; overrides = Dict("StorageDischargeDerate" => 1))
        @test vp(EP1, inputs, "A")[2] <= 5 + ATOL                 # 0.5 x 10 MW
        @test objective_value(EP1)≈0 atol=ATOL                    # B (30 MW, no derate) + A 5 >= 35
    end

    @testset "T5 refusals: Min_Duration != Max_Duration; inconsistent members" begin
        bad = ab_case(demand = [0.0, 35.0], g = [1.0, 0.0],
            a_kw = (; mind = 3, maxd = 5), b_kw = (; mind = 3, maxd = 5))
        @test occursin("Min_Duration 3", errmsg_ops(() -> ops_inputs(bad, ops_setup(bad))))
        ok = ab_case(demand = [0.0, 35.0], g = [1.0, 0.0])
        @test errmsg_ops(() -> ops_inputs(ok, ops_setup(ok))) == ""
        badeff = ab_case(demand = [0.0, 35.0], g = [1.0, 0.0], b_kw = (; eff = 0.9))
        @test occursin("inconsistent members", errmsg_ops(() -> ops_inputs(badeff, ops_setup(badeff))))
        # the same members in no group are accepted
        nogrp = ab_case(demand = [0.0, 35.0], g = [1.0, 0.0], b_kw = (; eff = 0.9),
            group = ["None", "None"])
        @test errmsg_ops(() -> ops_inputs(nogrp, ops_setup(nogrp))) == ""
    end

    @testset "T5 a group with LDS members is refused with Benders 1, accepted with 0" begin
        lds = ab_case(demand = [0.0, 35.0], g = [1.0, 0.0], a_kw = (; lds = 1), b_kw = (; lds = 1))
        @test occursin("not supported with Benders: 1",
            errmsg_ops(() -> ops_inputs(lds, ops_setup(lds; overrides = Dict("Benders" => 1)))))
        @test errmsg_ops(() -> ops_inputs(lds, ops_setup(lds; overrides = Dict("Benders" => 0)))) == ""
    end

    @testset "T5b #836 row with groups: substituted row, NSE 0" begin
        dir = ab_case(demand = [0.0, 0.0, 40.0, 40.0], g = [1.0, 1.0, 0.0, 0.0])
        EP, _ = ops_solve(dir; overrides = Dict("StorageChargeHeadroom" => 1))
        @test termination_status(EP) == MOI.OPTIMAL
        @test nse_mwh(EP)≈0 atol=ATOL               # the row left on eTotalCapEnergy[A]: NSE 40
        @test objective_value(EP)≈0 atol=ATOL
    end

    @testset "StorageChargeHeadroom 0 drops one row per storage and hour" begin
        dir = ab_case(demand = [0.0, 35.0], g = [1.0, 0.0], group = ["None", "None"])
        EP1, _ = ops_solve(dir; overrides = Dict("StorageChargeHeadroom" => 1))
        EP0, _ = ops_solve(dir; overrides = Dict("StorageChargeHeadroom" => 0))
        @test ncons(EP1) - ncons(EP0) == 2 * 2
        @test objective_value(EP0)≈objective_value(EP1) atol=ATOL
    end

    @testset "T10 a Stor_Op_Group column of None only: the model without the column" begin
        dir_none = ab_case(demand = [0.0, 35.0], g = [1.0, 0.0], group = ["None", "None"])
        dir_col = ops_case(; T = 2, demand = [0.0, 35.0], vre = [vre_row("G", 1, 100)],
            storage = [stor_row("A", 1, 10, 40), stor_row("B", 1, 30, 120)],
            variability = Dict("G" => [1.0, 0.0]))
        EPn, inn = ops_solve(dir_none)
        EPc, _ = ops_solve(dir_col)
        @test !inn["STOR_GROUPS_ACTIVE"]
        @test !has_name(EPn, :eTotalCapStorOp)
        @test ncons(EPn) == ncons(EPc)
        @test num_variables(EPn) == num_variables(EPc)
        @test objective_value(EPn) == objective_value(EPc)
    end
end

end # module
