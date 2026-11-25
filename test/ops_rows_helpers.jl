# Helpers for the tests of the operational rows ported from the chair's fork (reimagining-ra #204):
# test_ops_loaders.jl, test_ops_floors_budgets.jl, test_ops_storage_groups.jl, test_ops_hurdle.jl.
# `include` this file; it is not a test file.
#
# `ops_case` writes a tiny case to a fresh temporary folder. Conventions of the spec's tests
# (mwe/specs/base_ops_port.md section 9): omega_t = 1 (Sub_Weights = hours per period), one demand
# segment priced at VOLL = 1000 $/MWh, ParameterScale 0, UCommit 0, ObjScale 1, every resource
# existing (New_Build 0, Can_Retire 0), so the objective is operating cost only.
# Resource names must not parse as Booleans or numbers in CSV.jl ("T" is read as `true`).

using Test
using GenX
using JuMP, HiGHS
using Logging
using YAML

if isdefined(HiGHS, :Highs_resetGlobalScheduler)
    HiGHS.Highs_resetGlobalScheduler(1)
end

const VOLL = 1000.0

quiet_ops(f) = redirect_stdout(() -> with_logger(f, ConsoleLogger(stderr, Logging.Error)), devnull)

csvrow(xs) = join(string.(xs), ",")

const THERMAL_HDR = "Resource,Zone,Model,New_Build,Can_Retire,Existing_Cap_MW,Max_Cap_MW,Min_Cap_MW,Inv_Cost_per_MWyr,Fixed_OM_Cost_per_MWyr,Var_OM_Cost_per_MWh,Heat_Rate_MMBTU_per_MWh,Fuel,Cap_Size,Start_Cost_per_MW,Start_Fuel_MMBTU_per_MW,Up_Time,Down_Time,Ramp_Up_Percentage,Ramp_Dn_Percentage,Min_Power,Reg_Max,Rsv_Max,Reg_Cost,Rsv_Cost,region,cluster"
thermal_row(name, zone, cap, varom) =
    csvrow([name, zone, 2, 0, 0, cap, -1, 0, 0, 0, varom, 0, "None", 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, "z$zone", 1])

const VRE_HDR = "Resource,Zone,Num_VRE_Bins,New_Build,Can_Retire,Existing_Cap_MW,Max_Cap_MW,Min_Cap_MW,Inv_Cost_per_MWyr,Fixed_OM_Cost_per_MWyr,Var_OM_Cost_per_MWh,Reg_Max,Rsv_Max,Reg_Cost,Rsv_Cost,region,cluster"
vre_row(name, zone, cap) = csvrow([name, zone, 1, 0, 0, cap, -1, 0, 0, 0, 0, 0, 0, 0, 0, "z$zone", 1])

const HYDRO_HDR = "Resource,Zone,Hydro_Energy_To_Power_Ratio,Min_Power,Ramp_Up_Percentage,Ramp_Dn_Percentage,LDS,Existing_Cap_MW,Cap_Size,New_Build,Can_Retire,Min_Cap_MW,Max_Cap_MW,Inv_Cost_per_MWyr,Fixed_OM_Cost_per_MWyr,Var_OM_Cost_per_MWh,Eff_Up,Eff_Down,Reg_Max,Rsv_Max,Reg_Cost,Rsv_Cost,region,cluster"
hydro_row(name, zone, cap; lds = 0, min_power = 0) =
    csvrow([name, zone, 0, min_power, 1, 1, lds, cap, 1, 0, 0, 0, -1, 0, 0, 0, 1, 1, 0, 0, 0, 0, "z$zone", 1])

const STOR_HDR = "Resource,Zone,Model,LDS,New_Build,Can_Retire,Existing_Cap_MW,Existing_Cap_MWh,Max_Cap_MW,Max_Cap_MWh,Min_Cap_MW,Min_Cap_MWh,Inv_Cost_per_MWyr,Inv_Cost_per_MWhyr,Fixed_OM_Cost_per_MWyr,Fixed_OM_Cost_per_MWhyr,Var_OM_Cost_per_MWh,Var_OM_Cost_per_MWh_In,Self_Disch,Eff_Up,Eff_Down,Min_Duration,Max_Duration,Reg_Max,Rsv_Max,Reg_Cost,Rsv_Cost,region,cluster"
stor_row(name, zone, mw, mwh; lds = 0, model = 1, eff = 1, mind = mwh / mw, maxd = mwh / mw) =
    csvrow([name, zone, model, lds, 0, 0, mw, mwh, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, eff, eff,
        mind, maxd, 0, 0, 0, 0, "z$zone", 1])

"""
    ops_case(; T, demand, periods = 1, thermal = [], vre = [], hydro = [], storage = [],
             stor_group = nothing, variability = Dict(), pmin = nothing, budget = nothing,
             network = nothing, settings = Dict())

Write a case and return its folder. `demand`: Vector (1 zone) or Vector of Vectors (one per zone).
`thermal`/`vre`/`hydro`/`storage`: Vectors of resource rows (strings from the *_row helpers).
`stor_group`: Vector of group labels, one per storage row (adds a `Stor_Op_Group` column).
`variability`, `pmin`, `budget`: Dict(resource name => length-T Vector); `pmin`/`budget` =
`nothing` writes no file. `network`: (header, rows) of `Network.csv`. `settings` are merged into
genx_settings.yml.
"""
function ops_case(; T::Int, demand, periods::Int = 1, thermal = String[], vre = String[],
        hydro = String[], storage = String[], stor_group = nothing, variability = Dict(),
        pmin = nothing, budget = nothing, network = nothing, settings = Dict())
    dir = joinpath(mktempdir(), "case")
    for d in ("settings", "system", "resources", "policies")
        mkpath(joinpath(dir, d))
    end
    dem = demand isa Vector{<:Vector} ? demand : [demand]
    Z = length(dem)
    H = div(T, periods)
    open(joinpath(dir, "system", "Demand_data.csv"), "w") do io
        println(io, "Voll,Demand_Segment,Cost_of_Demand_Curtailment_per_MW,Max_Demand_Curtailment,",
            "Rep_Periods,Timesteps_per_Rep_Period,Sub_Weights,Time_Index,",
            join(["Demand_MW_z$z" for z in 1:Z], ","))
        for t in 1:T
            lead = t == 1 ? "$(Int(VOLL)),1,1,1,$periods,$H,$H" : ",,,,,,"
            lead = (t > 1 && t <= periods) ? ",,,,,,$H" : lead
            println(io, lead, ",", t, ",", join([dem[z][t] for z in 1:Z], ","))
        end
    end
    open(joinpath(dir, "system", "Fuels_data.csv"), "w") do io
        println(io, "Time_Index,None")
        println(io, "0,0.0")
        for t in 1:T
            println(io, t, ",0.0")
        end
    end
    names_all = String[]
    for (file, hdr, rows) in (("Thermal.csv", THERMAL_HDR, thermal), ("Vre.csv", VRE_HDR, vre),
        ("Hydro.csv", HYDRO_HDR, hydro), ("Storage.csv", STOR_HDR, storage))
        isempty(rows) && continue
        lines = copy(rows)
        h = hdr
        if file == "Storage.csv" && stor_group !== nothing
            h *= ",Stor_Op_Group"
            lines = [r * "," * string(stor_group[i]) for (i, r) in enumerate(rows)]
        end
        write(joinpath(dir, "resources", file), h * "\n" * join(lines, "\n") * "\n")
        append!(names_all, [first(split(r, ",")) for r in rows])
    end
    function write_profile(fname, prof, default)
        cols = sort(collect(keys(prof)))
        open(joinpath(dir, "system", fname), "w") do io
            println(io, join(["Time_Index"; cols], ","))
            for t in 1:T
                println(io, join([string(t); [string(prof[c][t]) for c in cols]], ","))
            end
        end
    end
    var = Dict{String, Vector{Float64}}()
    for n in names_all
        var[n] = haskey(variability, n) ? Float64.(variability[n]) : ones(T)
    end
    write_profile("Generators_variability.csv", var, 1.0)
    pmin === nothing || write_profile("Generators_Pmin.csv", pmin, 0.0)
    budget === nothing || write_profile("Hourly_energy_budget.csv", budget, 0.0)
    if network !== nothing
        hdr, rows = network
        write(joinpath(dir, "system", "Network.csv"), hdr * "\n" * join(rows, "\n") * "\n")
    end
    s = Dict{String, Any}("ParameterScale" => 0, "UCommit" => 0, "ObjScale" => 1,
        "NetworkExpansion" => 0, "Trans_Loss_Segments" => 1, "CO2Cap" => 0,
        "StorageLosses" => 1, "WriteShadowPrices" => 0, "OverwriteResults" => 1)
    merge!(s, Dict{String, Any}(string(k) => v for (k, v) in settings))
    open(joinpath(dir, "settings", "genx_settings.yml"), "w") do io
        for (k, v) in s
            println(io, k, ": ", v)
        end
    end
    write(joinpath(dir, "settings", "highs_settings.yml"),
        "Feasib_Tol: 1.0e-9\nOptimal_Tol: 1.0e-9\nPre_Solve: choose\nMethod: simplex\n")
    return dir
end

"Settings of the case folder, merged into the defaults, plus `overrides`."
function ops_setup(dir; overrides = Dict())
    setup = GenX.default_settings()
    merge!(setup, YAML.load_file(joinpath(dir, "settings", "genx_settings.yml")))
    merge!(setup, Dict{Any, Any}(string(k) => v for (k, v) in overrides))
    setup["EnableJuMPStringNames"] = true
    return setup
end

ops_inputs(dir, setup) = quiet_ops(() -> GenX.load_inputs(setup, dir))

"Load, build and solve; return (EP, inputs). Does not assert optimality."
function ops_solve(dir; overrides = Dict())
    setup = ops_setup(dir; overrides = overrides)
    inputs = ops_inputs(dir, setup)
    EP = quiet_ops() do
        OPT = GenX.configure_solver(joinpath(dir, "settings"), HiGHS.Optimizer)
        EP = GenX.generate_model(setup, inputs, OPT)
        set_silent(EP)
        set_attribute(EP, "threads", 1)
        optimize!(EP)
        EP
    end
    return EP, inputs
end

"Total non-served energy, MWh (omega = 1 in these cases, so unweighted)."
nse_mwh(EP) = sum(value.(EP[:vNSE]))

"Output of resource `name` over all hours."
function vp(EP, inputs, name)
    g = findfirst(==(name), inputs["RESOURCE_NAMES"])
    return [value(EP[:vP][g, t]) for t in 1:inputs["T"]]
end

"The message of the error thrown by `f()`, or \"\" if it throws none."
errmsg_ops(f) =
    try
        f()
        ""
    catch e
        sprint(showerror, e)
    end

has_name(EP, s::Symbol) = haskey(object_dictionary(EP), s)
