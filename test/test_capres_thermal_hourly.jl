# test/test_capres_thermal_hourly.jl — the thermal term of the capacity reserve margin row at hourly
# availability (setting CapResThermalHourlyAvailability) and the reserve-margin writers under storage
# operational groups (#204). Self-contained (GenX, JuMP, HiGHS, CSV, DataFrames, Test only): run it as
#   julia --project=<GenX clone> test/test_capres_thermal_hourly.jl
# Toy case: a copy of test/three_zones in a temporary folder, with
#   - MA_natural_gas_combined_cycle at availability 0.8 in hours 1-24 of Generators_variability.csv (1.0 before);
#   - a second battery MA_battery_2 (a copy of MA_battery), both in storage operational group MA_grp;
#   - one reserve-margin requirement of 0.15 over all three zones, derates thermal 0.93, VRE 0.8, storage 0.95.
module TestCapResThermalHourly

using Test
using GenX, JuMP, HiGHS, CSV, DataFrames

const SRC = joinpath(@__DIR__, "three_zones")
const DERATE = Dict("thermal" => 0.93, "vre" => 0.8, "storage" => 0.95)
const DOWN_HOURS = 1:24
const DOWN_RHO = 0.8

function make_case()
    case = joinpath(mktempdir(), "case")
    cp(SRC, case)
    gv_path = joinpath(case, "system", "Generators_variability.csv")
    gv = CSV.read(gv_path, DataFrame)
    gv[DOWN_HOURS, :MA_natural_gas_combined_cycle] .= DOWN_RHO
    gv[!, :MA_battery_2] = copy(gv[!, :MA_battery])
    CSV.write(gv_path, gv)

    st_path = joinpath(case, "resources", "Storage.csv")
    st = CSV.read(st_path, DataFrame; stringtype = String)
    @assert st[1, :Resource] == "MA_battery"
    push!(st, merge(copy(st[1, :]), (Resource = "MA_battery_2",)))
    st[!, :Stor_Op_Group] = [r == "MA_battery" || r == "MA_battery_2" ? "MA_grp" : "None"
                             for r in st.Resource]
    CSV.write(st_path, st)

    open(joinpath(case, "policies", "Capacity_reserve_margin.csv"), "w") do io
        write(io, ",Network_zones,CapRes_1\nMA,z1,0.15\nCT,z2,0.15\nME,z3,0.15\n")
    end
    lines = ["Resource,Derating_Factor_1"]
    for (file, cls) in (("Thermal.csv", "thermal"), ("Vre.csv", "vre"), ("Storage.csv", "storage"))
        df = CSV.read(joinpath(case, "resources", file), DataFrame; stringtype = String)
        append!(lines, ["$(r),$(DERATE[cls])" for r in df.Resource])
    end
    open(joinpath(case, "resources", "policy_assignments",
        "Resource_capacity_reserve_margin.csv"), "w") do io
        write(io, join(lines, "\n") * "\n")
    end
    return case
end

function toy_setup(extra::Dict = Dict())
    setup = GenX.default_settings()
    merge!(setup, Dict("NetworkExpansion" => 1, "Trans_Loss_Segments" => 1, "CO2Cap" => 0,
        "StorageLosses" => 1, "MinCapReq" => 1, "ParameterScale" => 0, "UCommit" => 0,
        "CapacityReserveMargin" => 1, "CapResELCC" => 0, "Hourly_Pmin" => 0,
        "Fixed_Dispatch" => 0, "WriteShadowPrices" => 1))
    merge!(setup, extra)
    setup["WriteOutputsSettingsDict"] = GenX.default_writeoutput()
    return setup
end

function build(case, setup)
    inputs = redirect_stdout(devnull) do
        GenX.load_inputs(setup, case)
    end
    OPT = GenX.configure_solver(case, HiGHS.Optimizer)
    EP = Model(OPT)
    redirect_stdout(devnull) do
        GenX.generate_model!(EP, setup, inputs)
    end
    return EP, inputs
end

# coefficient of the thermal resource's capacity in the row, per unit of eTotalCap
function thermal_row_coefficient(EP, y, t)
    v = EP[:vCAP][y]
    return coefficient(EP[:eCapResMarBalance][1, t], v) / coefficient(EP[:eTotalCap][y], v)
end

const CASE = make_case()

@testset "capres thermal hourly (#204)" begin
@testset "CapResThermalHourlyAvailability" begin
    EP1, inputs = build(CASE, toy_setup(Dict("CapResThermalHourlyAvailability" => 1)))
    gen = inputs["RESOURCES"]
    T = inputs["T"]
    y = findfirst(==("MA_natural_gas_combined_cycle"), inputs["RESOURCE_NAMES"])
    @test inputs["pP_Max"][y, 1] == DOWN_RHO          # the toy case has availability below 1
    @test inputs["pP_Max"][y, T] == 1.0
    @testset "key = 1: thermal coefficient in hour t = derate x pP_Max[y, t]" begin
        for t in (1, 24, 25, T)
            @test thermal_row_coefficient(EP1, y, t) ≈
                  GenX.derating_factor(gen[y], tag = 1) * inputs["pP_Max"][y, t]
        end
    end
    @testset "key = 0 and key absent: coefficient = derate in every hour, same row" begin
        EP0, _ = build(CASE, toy_setup(Dict("CapResThermalHourlyAvailability" => 0)))
        setup_absent = toy_setup()
        delete!(setup_absent, "CapResThermalHourlyAvailability")
        EPa, _ = build(CASE, setup_absent)
        for t in (1, 24, 25, T)
            @test thermal_row_coefficient(EP0, y, t) == GenX.derating_factor(gen[y], tag = 1)
            @test thermal_row_coefficient(EPa, y, t) == GenX.derating_factor(gen[y], tag = 1)
            @test string(EP0[:eCapResMarBalance][1, t]) == string(EPa[:eCapResMarBalance][1, t])
        end
    end
end

@testset "reserve-margin writers under storage operational groups" begin
    for key in (0, 1)
        setup = toy_setup(Dict("CapResThermalHourlyAvailability" => key))
        EP, inputs = build(CASE, setup)
        redirect_stdout(devnull) do
            GenX.solve_model(EP, setup)
        end
        @test termination_status(EP) == MOI.OPTIMAL
        @test length(inputs["STOR_OPERATIONAL"]) < length(inputs["STOR_ALL"])   # groups are on
        out = mktempdir()
        df = GenX.write_reserve_margin_revenue(out, inputs, setup, EP)
        @test isfile(joinpath(out, "ReserveMarginRevenue.csv"))
        # thermal revenue = sum_t derate x (pP_Max[y,t] if key 1 else 1) x capacity x dual_t
        y = findfirst(==("MA_natural_gas_combined_cycle"), inputs["RESOURCE_NAMES"])
        price = GenX.capacity_reserve_margin_price(EP, inputs, setup, 1) .* inputs["omega"]
        rho = key == 1 ? inputs["pP_Max"][y, :] : ones(inputs["T"])
        expected = GenX.derating_factor(inputs["RESOURCES"][y], tag = 1) *
                   value(EP[:eTotalCap][y]) * sum(rho .* price)
        @test df[y, :CapRes_1] ≈ expected rtol = 1e-9
        GenX.write_capacity_value(out, inputs, setup, EP)
        @test isfile(joinpath(out, "CapacityValue.csv"))
        res = redirect_stdout(devnull) do
            GenX.write_outputs(EP, joinpath(out, "full"), setup, inputs)
        end
        @test isfile(joinpath(res, "ReserveMarginRevenue.csv"))
        @test isfile(joinpath(res, "NetRevenue.csv"))
    end
end
end # outer testset

end # module TestCapResThermalHourly
