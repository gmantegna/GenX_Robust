@doc raw"""
    load_nse_disabled_periods!(setup::Dict, path::AbstractString, inputs::Dict)

Read `NSE_disabled_periods.csv` from the policies folder (`NSEDisabledPeriods: 1`): one column
`Rep_Period`, one row per representative-period index in which non-served energy is forbidden
(every hour of the period gets `vNSE[s,t,z] == 0`, see `non_served_energy!`). Stores the sorted,
de-duplicated indices as `inputs["NSE_DISABLED_PERIODS"]`. An index outside
`1:inputs["REP_PERIOD"]` (including 0) is an error.

`Rep_Period` is the index `1:REP_PERIOD` of the representative period as the model numbers it,
not the chronological period id of `Period_map.csv`.

With `Benders: 1` a listed period's subproblem loses complete recourse: with `vNSE` fixed to zero
it can be infeasible for a given investment plan, so Benders has to be able to return feasibility
cuts (`ExpectFeasibleSubproblems` must not be `true`).
"""
function load_nse_disabled_periods!(setup::Dict, path::AbstractString, inputs::Dict)
    filename = "NSE_disabled_periods.csv"
    df = load_dataframe(joinpath(path, filename))
    "Rep_Period" in names(df) ||
        error("$filename must have a column Rep_Period (got $(names(df))).")
    periods = sort(unique(Int.(df[!, :Rep_Period])))
    W = inputs["REP_PERIOD"]
    bad = filter(p -> !(1 <= p <= W), periods)
    isempty(bad) ||
        error("$filename lists representative period(s) $bad outside 1:$W.")
    inputs["NSE_DISABLED_PERIODS"] = periods
    println(filename * " Successfully Read!")
end

@doc raw"""
    nse_disabled_hours(inputs::Dict) -> Vector{Int}

Hours `t` of the model in which non-served energy is disabled. For the full (monolithic) inputs
these are the hours `hours_per_subperiod × (p − 1) + 1 … hours_per_subperiod × p` of every listed
period `p`; for a Benders subproblem's sliced inputs (`inputs["SubPeriod"]` set, `T` = one
period) they are all of its hours if its period is listed, else none. Empty when the setting is
off or the list is empty.
"""
function nse_disabled_hours(inputs::Dict)
    periods = get(inputs, "NSE_DISABLED_PERIODS", Int[])
    isempty(periods) && return Int[]
    if haskey(inputs, "SubPeriod")
        return inputs["SubPeriod"] in periods ? collect(1:inputs["T"]) : Int[]
    end
    H = inputs["hours_per_subperiod"]
    hours = Int[]
    for p in periods
        append!(hours, (H * (p - 1) + 1):(H * p))
    end
    return hours
end
