module TestPRMELCC

# Tests of the PRM-with-ELCC-surfaces policy (`PRM_ELCC`, src/model/policies/prm_elcc.jl).
# Case: test/prm_elcc/1_zone_prm (1 zone, 1 representative period of 24 h weighted 8760 so that
# omega_t = 365, constant demand D = 1000 MW, VoLL 50000 $/MWh in one segment, ParameterScale 1,
# plain LP). Three resources:
#   gas_ct      Thermal, new build, c_g = 80000 + 20000 = 100000 $/MW-yr, marginal cost
#               2 + 8 x 3.5 = 30 $/MWh; firm, NQC 0.9.
#   solar_pv    Vre, FIXED existing 300 MW (New_Build 0, Can_Retire 0), FOM 10000 $/MW-yr,
#               availability 0 at night and a daytime bell with sum_t a_t = 7.8, max 1.0; on the
#               surface axis `solar` (multiplier 1).
#   battery_4h  Storage, new build, Min_Duration = Max_Duration = 4, c_b = 40000 + 4000 +
#               4 x 4000 = 60000 $/MW-yr, round-trip efficiency 0.09; on the axis `storage`.
# Surface `sol_stor`, two facets:  f1: ELCC <= 0   + 0.5 solar + 1.00 storage
#                                  f2: ELCC <= 200 + 0.1 solar + 0.75 storage
# TRN = 1700 MW.
#
# Hand solution. Energy: at night solar is 0 and the battery is useless for energy (shaving 1 MW
# of gas capacity over the 12 night hours would burn 12/0.09 MWh of gas a day, 1.46 M$/MW-yr >>
# c_g; shedding costs VoLL), so x_g = D = 1000 MW and the only free decision is the battery
# built for the PRM. Firm accreditation 0.9 x 1000 = 900, so the surface must deliver
# E = 1700 - 900 = 800 MW. With x_s = 300 fixed: f1 = 150 + x_b, f2 = 230 + 0.75 x_b; they
# cross at x_b = 320 (the min switches inside the range), and E = 800 needs x_b = 570/0.75 =
# 760 MW on f2 (f1 is then 910 > 800: slack, f2 is the lower plane). Dual feasibility: the PRM
# multiplier is lambda = c_b / 0.75 = 80000 $/MW-yr; gas has reduced cost 100000 - 0.9 x 80000
# = 28000 > 0, so it is not built for the PRM. Annual cost:
#   gas fixed      100000 x 1000                      = 100,000,000
#   gas variable   30 x 365 x (24 x 1000 - 300 x 7.8) = 237,177,000
#   solar FOM      10000 x 300                        =   3,000,000
#   battery        60000 x 760                        =  45,600,000
#   total                                              385,777,000 $/yr  (385.777 M$ in model units)
# Facet f1 alone (f2 removed): x_b = 650 -> 379,177,000. No binding PRM (TRN tiny): x_b = 0 ->
# 340,177,000, dual 0.

using Test
using GenX
using JuMP, HiGHS
using Logging
using YAML, CSV, DataFrames

const CASE = joinpath(@__DIR__, "prm_elcc", "1_zone_prm")

const D_MW = 1000.0
const TRN_MW = 1700.0
const XB_HAND = 760.0
const OBJ_HAND = 385_777_000.0          # $/yr
const OBJ_F1_ONLY = 379_177_000.0       # $/yr, facet f2 removed
const OBJ_NO_PRM = 340_177_000.0        # $/yr, TRN tiny
const LAMBDA_HAND = 80_000.0            # $/MW-yr
const RTOL = 1e-6

quiet(f) = redirect_stdout(() -> with_logger(f, ConsoleLogger(stderr, Logging.Warn)), devnull)

"""
Copy of the case in a temporary folder, with settings `overrides` merged into genx_settings.yml,
`trn` written to policies/PRM_ELCC.csv and the facets restricted to `facets` (nothing = all).
The repo's case folder is never written to.
"""
function scratch_case(; overrides::Dict = Dict(), trn = TRN_MW, facets = nothing)
    dir = joinpath(mktempdir(), "case")
    cp(CASE, dir)
    for d in filter(startswith("results"), readdir(dir))
        rm(joinpath(dir, d); recursive = true)
    end
    f = joinpath(dir, "settings", "genx_settings.yml")
    y = YAML.load_file(f)
    merge!(y, overrides)
    YAML.write_file(f, y)
    write(joinpath(dir, "policies", "PRM_ELCC.csv"), "TRN_MW\n$trn\n")
    if facets !== nothing
        df = CSV.read(joinpath(dir, "policies", "PRM_ELCC_facets.csv"), DataFrame)
        CSV.write(joinpath(dir, "policies", "PRM_ELCC_facets.csv"), df[in.(df.Facet, Ref(facets)), :])
    end
    return dir
end

scale(dir) = get(YAML.load_file(joinpath(dir, "settings", "genx_settings.yml")), "ParameterScale", 0) == 1 ?
             GenX.ModelScalingFactor : 1.0

"Run the case end to end; return (objective in \$/yr, capacity DataFrame, results dir)."
function run_case(dir)
    quiet(() -> GenX.run_genx_case!(dir, HiGHS.Optimizer))
    out = joinpath(dir, "results")
    status = CSV.read(joinpath(out, "status.csv"), DataFrame)
    cap = CSV.read(joinpath(out, "capacity.csv"), DataFrame)
    obj = status.Objval[1] * scale(dir)^2       # M\$ -> \$ when ParameterScale = 1
    return obj, cap, out
end

endcap(cap, r) = cap[cap.Resource .== r, :EndCap][1]

function read_prm(out)
    f = joinpath(out, "prm_elcc.csv")
    @test isfile(f)
    return isfile(f) ? CSV.read(f, DataFrame; missingstring = "") : DataFrame()
end

item(df, it; surface = missing, facet = missing, resource = missing, axis = missing) = begin
    m = df.Item .== it
    ismissing(surface) || (m .&= coalesce.(df.Surface .== surface, false))
    ismissing(facet) || (m .&= coalesce.(df.Facet .== facet, false))
    ismissing(resource) || (m .&= coalesce.(df.Resource .== resource, false))
    ismissing(axis) || (m .&= coalesce.(df.Axis .== axis, false))
    df[m, :]
end

@testset "PRM_ELCC" begin
    @testset "T0 setting off: model unchanged, defaults" begin
        d = GenX.default_settings()
        @test haskey(d, "PRM_ELCC") && d["PRM_ELCC"] == 0
        dir = scratch_case(overrides = Dict("PRM_ELCC" => 0))
        obj, cap, out = run_case(dir)
        @test obj≈OBJ_NO_PRM rtol=RTOL
        @test endcap(cap, "battery_4h")≈0.0 atol=1e-6
        @test !isfile(joinpath(out, "prm_elcc.csv"))
    end

    @testset "T1 optimum matches the hand solution (ParameterScale 1 and 0)" begin
        for ps in (1, 0)
            dir = scratch_case(overrides = Dict("ParameterScale" => ps))
            obj, cap, out = run_case(dir)
            println("T1 ParameterScale=$ps: objective = $obj \$/yr (hand $OBJ_HAND), " *
                    "gas = $(endcap(cap, "gas_ct")), battery = $(endcap(cap, "battery_4h")), " *
                    "solar = $(endcap(cap, "solar_pv")) MW")
            @test obj≈OBJ_HAND rtol=RTOL
            @test endcap(cap, "gas_ct")≈D_MW rtol=RTOL
            @test endcap(cap, "battery_4h")≈XB_HAND rtol=RTOL
            @test endcap(cap, "solar_pv")≈300.0 rtol=RTOL
            # regime assumptions of the hand solution: no shed, no curtailment, battery idle
            nse = CSV.read(joinpath(out, "nse.csv"), DataFrame; header = false)
            @test nse[3, 1] == "AnnualSum"
            @test abs(parse(Float64, string(nse[3, end]))) < 1e-6
            curt = CSV.read(joinpath(out, "curtail.csv"), DataFrame; header = false)
            @test all(abs(parse(Float64, string(x))) < 1e-6 for x in curt[3, 2:end])
            ch = CSV.read(joinpath(out, "charge.csv"), DataFrame; header = false)
            @test all(abs(parse(Float64, string(x))) < 1e-6 for x in ch[3, 2:end])

            prm = read_prm(out)
            if nrow(prm) > 0
                # dual of cPRM_ELCC in $/MW-yr, independent of ParameterScale
                @test item(prm, "Dual").Value[1]≈LAMBDA_HAND rtol=RTOL
                @test item(prm, "TRN").Value[1]≈TRN_MW rtol=RTOL
                @test item(prm, "SurfaceELCC"; surface = "sol_stor").Value[1]≈TRN_MW - 0.9 * D_MW rtol=RTOL
                @test item(prm, "FirmAccredited"; resource = "gas_ct").Value[1]≈0.9 * D_MW rtol=RTOL
                @test item(prm, "TotalAccredited").Value[1]≈TRN_MW rtol=RTOL
                @test item(prm, "AxisCapacity"; axis = "storage").Value[1]≈XB_HAND rtol=RTOL
                @test item(prm, "AxisCapacity"; axis = "solar").Value[1]≈300.0 rtol=RTOL
            end
        end
    end

    @testset "T2 the binding facet is the lower plane" begin
        dir = scratch_case()
        obj, cap, out = run_case(dir)
        prm = read_prm(out)
        if nrow(prm) > 0
            b = item(prm, "BindingFacet"; surface = "sol_stor")
            @test nrow(b) == 1
            @test b.Facet[1] == "f2"
            @test abs(b.Value[1]) < 1e-6                     # zero slack on f2
            @test item(prm, "FacetSlack"; facet = "f1").Value[1]≈110.0 rtol=1e-5   # 910 - 800
            @test item(prm, "FacetDual"; facet = "f2").Value[1]≈LAMBDA_HAND rtol=RTOL
            @test abs(item(prm, "FacetDual"; facet = "f1").Value[1]) < 1e-6
            @test item(prm, "FacetValue"; facet = "f2").Value[1]≈800.0 rtol=RTOL
            @test item(prm, "FacetValue"; facet = "f1").Value[1]≈910.0 rtol=RTOL
        end
    end

    @testset "T3 facets bind: removing the binding facet changes the optimum" begin
        obj_f1, cap_f1, _ = run_case(scratch_case(facets = ["f1"]))
        println("T3 f1 only: objective = $obj_f1 (hand $OBJ_F1_ONLY), battery = $(endcap(cap_f1, "battery_4h"))")
        @test obj_f1≈OBJ_F1_ONLY rtol=RTOL
        @test endcap(cap_f1, "battery_4h")≈650.0 rtol=RTOL
        @test !(obj_f1 ≈ OBJ_HAND)
        # the slack facet can go without changing anything
        obj_f2, cap_f2, _ = run_case(scratch_case(facets = ["f2"]))
        @test obj_f2≈OBJ_HAND rtol=RTOL
        @test endcap(cap_f2, "battery_4h")≈XB_HAND rtol=RTOL
    end

    @testset "T4 prm_elcc.csv: positive dual when binding, zero when TRN is tiny" begin
        _, _, out = run_case(scratch_case())
        prm = read_prm(out)
        nrow(prm) > 0 && @test item(prm, "Dual").Value[1] > 0
        obj0, cap0, out0 = run_case(scratch_case(trn = 1.0))
        @test obj0≈OBJ_NO_PRM rtol=RTOL
        @test endcap(cap0, "battery_4h")≈0.0 atol=1e-6
        prm0 = read_prm(out0)
        if nrow(prm0) > 0
            @test abs(item(prm0, "Dual").Value[1]) < 1e-9
            @test item(prm0, "TotalAccredited").Value[1] >= 1.0 - 1e-6
        end
    end

    @testset "T5 model structure and input validation" begin
        setup = quiet() do
            GenX.configure_settings(joinpath(CASE, "settings", "genx_settings.yml"),
                joinpath(CASE, "settings", "output_settings.yml"))
        end
        inputs = quiet(() -> GenX.load_inputs(setup, CASE))
        EP = quiet(() -> GenX.generate_model(setup, inputs,
            GenX.configure_solver(joinpath(CASE, "settings"), HiGHS.Optimizer)))
        od = object_dictionary(EP)
        @test haskey(od, :cPRM_ELCC) && haskey(od, :cELCCFacet) && haskey(od, :vELCC)
        @test length(EP[:vELCC]) == 1
        @test !has_lower_bound(EP[:vELCC]["sol_stor"]) && !has_upper_bound(EP[:vELCC]["sol_stor"])
        @test constraint_object(EP[:cPRM_ELCC]).set isa MOI.GreaterThan
        @test normalized_rhs(EP[:cPRM_ELCC]) ≈ TRN_MW / GenX.ModelScalingFactor
        facets = collect(keys(EP[:cELCCFacet].data))
        @test length(facets) == 2
        for k in facets
            c = EP[:cELCCFacet][k...]
            @test constraint_object(c).set isa MOI.LessThan
            @test normalized_coefficient(c, EP[:vELCC]["sol_stor"]) ≈ 1.0
        end
        c2 = EP[:cELCCFacet]["sol_stor", "f2"]
        @test normalized_rhs(c2) ≈ 200.0 / GenX.ModelScalingFactor + 0.1 * 300.0 / GenX.ModelScalingFactor
        # storage counts by POWER capacity: the energy variable has no coefficient
        @test normalized_coefficient(c2, EP[:vCAPENERGY][inputs["STOR_ALL"][1]]) == 0.0
        # Benders is refused, with a clear message
        sb = copy(setup)
        sb["Benders"] = 1
        err = try
            GenX.generate_model(sb, inputs, HiGHS.Optimizer)
            "no error"
        catch e
            e isa ErrorException ? e.msg : "not an ErrorException: $(typeof(e))"
        end
        @test occursin("PRM_ELCC is a planning-level constraint; Benders support not implemented", err)
        # a resource missing from Resource_prm_elcc.csv contributes nothing, with a warning
        dir = scratch_case()
        rf = joinpath(dir, "resources", "policy_assignments", "Resource_prm_elcc.csv")
        write(rf, "Resource,PRM_ELCC_NQC,PRM_ELCC_Axis,PRM_ELCC_Axis_Multiplier\n" *
                  "solar_pv,,solar,1.0\nbattery_4h,,storage,1.0\n")
        s2 = quiet(() -> GenX.configure_settings(joinpath(dir, "settings", "genx_settings.yml"),
            joinpath(dir, "settings", "output_settings.yml")))
        @test_logs (:warn, r"gas_ct") match_mode=:any redirect_stdout(() -> GenX.load_inputs(s2, dir), devnull)
    end
end

end # module TestPRMELCC
