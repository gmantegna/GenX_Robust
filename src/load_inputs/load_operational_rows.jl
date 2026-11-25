@doc raw"""
    load_hourly_pmin!(setup::Dict, path::AbstractString, inputs::Dict)

With `Hourly_Pmin: 1`, reads `Generators_Pmin.csv` (same folder as `Generators_variability.csv`;
the file is required): one column per resource, the hourly minimum output as a fraction of total
capacity. Stores `inputs["pP_Min"]` (G x T; a resource without a column gets 0) and
`inputs["PMIN_RIDS"]` (the resources with a column, the only ones that get a floor row).
Ported from the chair's GenX fork (`load_generators_pmin!`).
"""
function load_hourly_pmin!(setup::Dict, path::AbstractString, inputs::Dict)
    get(setup, "Hourly_Pmin", 0) == 1 || return nothing
    my_dir = get_systemfiles_path(setup, joinpath(path, setup["TimeDomainReductionFolder"]), path)
    filename = joinpath(my_dir, "Generators_Pmin.csv")
    isfile(filename) || error("Hourly_Pmin is 1 but $filename does not exist.")
    inputs["pP_Min"], inputs["PMIN_RIDS"] = _resource_profile_matrix(load_dataframe(filename),
        inputs, "Generators_Pmin.csv")
    println("Generators_Pmin.csv Successfully Read!")
    return nothing
end

@doc raw"""
    load_energy_budget!(setup::Dict, path::AbstractString, inputs::Dict)

If `Hourly_energy_budget.csv` exists in the system folder (the chair's fork's switch), reads it:
one column per budgeted resource, the hourly energy budget as a fraction of total capacity
(MWh per MW per hour). Stores `inputs["pBudget"]` (G x T; 0 without a column) and
`inputs["BUDGET_RIDS"]` (the resources with a column).
"""
function load_energy_budget!(setup::Dict, path::AbstractString, inputs::Dict)
    isfile(joinpath(path, setup["SystemFolder"], "Hourly_energy_budget.csv")) || return nothing
    my_dir = get_systemfiles_path(setup, joinpath(path, setup["TimeDomainReductionFolder"]), path)
    filename = joinpath(my_dir, "Hourly_energy_budget.csv")
    inputs["pBudget"], inputs["BUDGET_RIDS"] = _resource_profile_matrix(load_dataframe(filename),
        inputs, "Hourly_energy_budget.csv")
    println("Hourly_energy_budget.csv Successfully Read!")
    return nothing
end

function _resource_profile_matrix(df::DataFrame, inputs::Dict, filename::AbstractString)
    G, T = inputs["G"], inputs["T"]
    nrow(df) >= T || error("$filename has $(nrow(df)) rows; at least T = $T are needed.")
    names_df = names(df)
    M = zeros(Float64, G, T)
    rids = Int[]
    for (g, r) in enumerate(inputs["RESOURCE_NAMES"])
        if r in names_df
            M[g, :] = Float64.(df[1:T, r])
            push!(rids, g)
        end
    end
    return M, rids
end

@doc raw"""
    operational_rows_voll(inputs::Dict)

The largest non-served-energy segment cost over segments (`Voll` x
`Cost_of_Demand_Curtailment_per_MW`), the VoLL that `OperationalFloorPenalty` must exceed.
"""
operational_rows_voll(inputs::Dict) = maximum(inputs["pC_D_Curtail"])

@doc raw"""
    validate_operational_rows!(setup::Dict, inputs::Dict)

Refusals for the operational rows ported from the chair's fork (reimagining-ra #204):
`Fixed_Dispatch: 1` (not ported); `ParameterScale: 1` with any of these inputs or keys in use
(the rows are not scaled); `OperationalFloorPenalty` not a number above VoLL; an
`EnergyBudgetWindowHours` that is not a positive integer dividing the subperiod length; a
budgeted hydro resource with `LDS: 1`; a 0/1 key with another value.
"""
function validate_operational_rows!(setup::Dict, inputs::Dict)
    for k in ("Hourly_Pmin", "Fixed_Dispatch", "HydroBudgetMaxProfile",
        "StorageDischargeDerate", "StorageChargeHeadroom")
        get(setup, k, 0) in (0, 1) || error("$k must be 0 or 1 (got $(setup[k])).")
    end
    get(setup, "Fixed_Dispatch", 0) == 1 &&
        error("Fixed_Dispatch: 1 (Generators_fixed_dispatch.csv) is not supported on this branch.")

    in_use = haskey(inputs, "pP_Min") || haskey(inputs, "pBudget") ||
             haskey(setup, "OperationalFloorPenalty") ||
             haskey(setup, "EnergyBudgetWindowHours") ||
             get(setup, "HydroBudgetMaxProfile", 0) == 1 ||
             get(setup, "StorageDischargeDerate", 0) == 1
    if in_use && setup["ParameterScale"] == 1
        error("ParameterScale: 1 is not supported with hourly Pmin, energy budgets, " *
              "OperationalFloorPenalty, EnergyBudgetWindowHours, HydroBudgetMaxProfile or " *
              "StorageDischargeDerate (these rows are not scaled).")
    end

    if haskey(setup, "OperationalFloorPenalty")
        P = setup["OperationalFloorPenalty"]
        voll = operational_rows_voll(inputs)
        (P isa Real && isfinite(P) && P > voll) ||
            error("OperationalFloorPenalty must be a number above VoLL = $voll (the largest " *
                  "non-served-energy segment cost); got $P.")
    end

    H = inputs["hours_per_subperiod"]
    if haskey(setup, "EnergyBudgetWindowHours")
        W = setup["EnergyBudgetWindowHours"]
        (W isa Integer && W >= 1 && H % W == 0) ||
            error("EnergyBudgetWindowHours must be a positive integer that divides the " *
                  "subperiod length $H; got $W.")
    end

    if haskey(inputs, "BUDGET_RIDS")
        bad = intersect(inputs["BUDGET_RIDS"], inputs["STOR_HYDRO_LONG_DURATION"])
        isempty(bad) ||
            error("Hydro resources with an energy budget column cannot have LDS = 1: " *
                  join(inputs["RESOURCE_NAMES"][bad], ", "))
    end
    return nothing
end

"The energy budget window in hours: `EnergyBudgetWindowHours`, default the subperiod length."
energy_budget_window(setup::Dict, inputs::Dict) =
    get(setup, "EnergyBudgetWindowHours", inputs["hours_per_subperiod"])
