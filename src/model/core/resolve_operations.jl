@doc raw"""
    resolve_operations!(EP::Model, inputs::Dict, setup::Dict)

Operational rows ported from the chair's GenX fork (reimagining-ra #204). Built in
`operation_model!`, so in the monolithic model and in every Benders subproblem. With neither
`Generators_Pmin.csv` in use nor `Hourly_energy_budget.csv` present, nothing is built.

**Hourly minimum output** (`Hourly_Pmin: 1`), for every resource $g$ with a column in
`Generators_Pmin.csv` and every hour $t$ (`cHourlyPmin`):
```math
\Theta_{g,t} \;(+\, u_{g,t}) \geq \pi_{g,t}\,\Delta_g
```
**Energy budget** (file `Hourly_energy_budget.csv`), for every resource $g$ with a column and every
window $w = [s, s+W-1]$, $W$ = `EnergyBudgetWindowHours` (default: the subperiod length), windows
starting at each subperiod start and every $W$ hours after it (`cHourlyBudget`):
```math
\sum_{\tau \in w} \Theta_{g,\tau} \;(+\, v_{g,w}) = \Delta_g \sum_{\tau \in w} e_{g,\tau}
```
The slacks $u \geq 0$ (`vFloorSlack`) and $v \geq 0$ (`vBudgetSpill`) exist only when
`OperationalFloorPenalty` $P$ is set; each MWh costs $P\,\omega_t$ ($\omega$ of the hour, or of the
window's first hour), in `eTotalCFloorSlack` and `eTotalCBudgetSpill`, added to the objective.
Without the key the rows are hard. There is no upward budget slack.

**Hourly hydro maximum** (`HydroBudgetMaxProfile: 1`), for budgeted hydro only
(`cHydroBudgetMaxProfile`): $\Theta_{g,t} \leq \rho_{g,t}\,\Delta_g$, $\rho$ = the resource's
`Generators_variability.csv` column. The reservoir rows of budgeted hydro, which read that column
as inflow, are not built (`hydro_res!`).
"""
function resolve_operations!(EP::Model, inputs::Dict, setup::Dict)
    has_pmin = haskey(inputs, "pP_Min")
    has_budget = haskey(inputs, "pBudget")
    (has_pmin || has_budget) || return nothing
    println("Hourly Pmin and Energy Budget Module")

    T = inputs["T"]
    omega = inputs["omega"]
    vP = EP[:vP]
    eTotalCap = EP[:eTotalCap]
    slack = haskey(setup, "OperationalFloorPenalty")
    P = slack ? Float64(setup["OperationalFloorPenalty"]) : 0.0

    if has_pmin
        PMIN = inputs["PMIN_RIDS"]
        pmin = inputs["pP_Min"]
        if slack
            @variable(EP, vFloorSlack[g in PMIN, t = 1:T]>=0)
            @constraint(EP, cHourlyPmin[g in PMIN, t = 1:T],
                vP[g, t] + vFloorSlack[g, t]>=pmin[g, t] * eTotalCap[g])
            @expression(EP, eTotalCFloorSlack,
                sum(P * omega[t] * vFloorSlack[g, t] for g in PMIN, t in 1:T))
            add_to_expression!(EP[:eObj], eTotalCFloorSlack)
        else
            @constraint(EP, cHourlyPmin[g in PMIN, t = 1:T],
                vP[g, t]>=pmin[g, t] * eTotalCap[g])
        end
    end

    if has_budget
        BUDGET = inputs["BUDGET_RIDS"]
        e = inputs["pBudget"]
        H = inputs["hours_per_subperiod"]
        W = energy_budget_window(setup, inputs)
        WSTART = [s + k * W for s in inputs["START_SUBPERIODS"] for k in 0:(div(H, W) - 1)]
        if slack
            @variable(EP, vBudgetSpill[g in BUDGET, w in WSTART]>=0)
            @constraint(EP, cHourlyBudget[w in WSTART, g in BUDGET],
                sum(vP[g, tau] for tau in w:(w + W - 1)) + vBudgetSpill[g, w]==
                eTotalCap[g] * sum(e[g, tau] for tau in w:(w + W - 1)))
            @expression(EP, eTotalCBudgetSpill,
                sum(P * omega[w] * vBudgetSpill[g, w] for g in BUDGET, w in WSTART))
            add_to_expression!(EP[:eObj], eTotalCBudgetSpill)
        else
            @constraint(EP, cHourlyBudget[w in WSTART, g in BUDGET],
                sum(vP[g, tau] for tau in w:(w + W - 1))==
                eTotalCap[g] * sum(e[g, tau] for tau in w:(w + W - 1)))
        end
        if get(setup, "HydroBudgetMaxProfile", 0) == 1
            BH = intersect(BUDGET, inputs["HYDRO_RES"])
            @constraint(EP, cHydroBudgetMaxProfile[g in BH, t = 1:T],
                vP[g, t]<=inputs["pP_Max"][g, t] * eTotalCap[g])
        end
    end
    return nothing
end

"Hydro resources with an energy budget column (their reservoir rows are not built)."
budgeted_hydro(inputs::Dict) =
    haskey(inputs, "BUDGET_RIDS") ? intersect(inputs["HYDRO_RES"], inputs["BUDGET_RIDS"]) : Int[]
