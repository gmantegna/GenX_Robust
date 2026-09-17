@doc raw"""
	nse_budget!(EP::Model, inputs::Dict, setup::Dict)

Hard cap on expected annual non-served energy (expected unserved energy, EUE), system-wide.
Activated by `NSEBudget: 1` in `genx_settings.yml`; the cap $\bar{E}$ is `NSEBudgetTargetMWh`
(MWh per year; divided by `ModelScalingFactor` internally when `ParameterScale: 1`).

**Monolithic model** (this function). One row over all segments, hours and zones:

```math
\begin{aligned}
	\sum_{z \in \mathcal{Z}} \sum_{t \in \mathcal{T}} \sum_{s \in \mathcal{S}} \omega_{t} \Lambda_{s,t,z} \leq \bar{E}
\end{aligned}
```

**Benders decomposition** ([`nse_budget_planning!`](@ref), [`nse_budget_subperiod!`](@ref)).
The device is the one of the CO$_2$ cap (`co2_cap_planning!` / `co2_cap_subperiod!`): the
planning problem carries one budget per representative period $w$ and the cap on their sum,

```math
\begin{aligned}
	q_w \geq 0 \quad \forall w \in \mathcal{W}, \qquad \sum_{w \in \mathcal{W}} q_w \leq \bar{E},
\end{aligned}
```

and the operational subproblem of period $w$ is held to its budget,

```math
\begin{aligned}
	\sum_{z \in \mathcal{Z}} \sum_{t \in \mathcal{T}_w} \sum_{s \in \mathcal{S}} \omega_{t} \Lambda_{s,t,z} \leq q_w .
\end{aligned}
```

$q_w$ is the variable `vNSEbudget[w]`. It exists under the same JuMP name in the planning
problem and in subproblem $w$, so it is a linking variable: `MacroEnergySolvers` fixes it in the
subproblem, and the dual of that fixing enters the ordinary optimality cut, which is therefore a
cut in (capacity, budget) space. A subproblem whose budget cannot be met at the trial capacities
is infeasible; with `ExpectFeasibleSubproblems: false` in `benders_settings.yml` it returns a
feasibility cut through the `slack_max` phase-I problem of `MacroEnergySolvers`, which relaxes
every affine row, this one included. `ExpectFeasibleSubproblems: true` is an error in this mode
because complete recourse does not hold.

**Weighting convention.** Every quantity here is in GenX's weighted convention:
$\omega_t = $ `Sub_Weights[w]` / `Timesteps_per_Rep_Period` for $t \in \mathcal{T}_w$, and the
subproblem objective, hence `vTHETA[w]`, carries the same $\omega_t$. When the `Sub_Weights` sum
to 8760 and period $w$ is a scenario-year of probability $\pi_w$ (8760 hourly steps,
`Sub_Weights[w]` $= 8760\,\pi_w$), $\omega_t = \pi_w$, the left-hand side of the subproblem row
is $\pi_w \times$ (annual shed of year $w$), `vNSEbudget[w]` is $\pi_w \times$ (the annual shed
budget of year $w$), and $\sum_w$ `vNSEbudget[w]` is expected annual unserved energy in MWh. The
cap is therefore a plain sum, with **no further probability weight**. Because objective and row
both carry $\omega_t$, the dual of the subproblem row and the dual of the cap are prices per
unit of energy (\$/MWh times `ObjScale`), not probability-weighted prices. If the `Sub_Weights`
do not sum to 8760 the cap is on $\sum_t \omega_t \Lambda_t$ whatever that sum represents; a
warning is issued.

**VoLL.** With `NSEBudgetRemoveVoLL: 1` (default) the cost of non-served energy is removed from
the objective in every segment, so the multiplier of the cap is the whole implied price of
unserved energy. If the VoLL form (single segment, price VoLL) has optimal value $z_{VoLL}$ with
expected shed $\bar{E}^*$, the budget form with $\bar{E} = \bar{E}^*$ has optimal value
$z_{VoLL} - \mathrm{VoLL}\,\bar{E}^*$. The identity is one of optimal values, not of portfolios
or returned duals, and it needs a single segment: with several segments the cap does not
distinguish them once their prices are removed (a warning is issued). With
`NSEBudgetRemoveVoLL: 0` the segment prices stay in the objective and the cap's multiplier is
an increment on top of them.
"""
function nse_budget!(EP::Model, inputs::Dict, setup::Dict)
    println("NSE Budget Policy Module")

    target = nse_budget_target(inputs, setup)

    # Expected NSE over the whole horizon, weighted convention (MWh, or GWh if ParameterScale = 1)
    @expression(EP, eNSEBudgetExpectedNSE, nse_budget_expected_nse(EP, inputs))

    @constraint(EP, cNSEBudget, eNSEBudgetExpectedNSE<=target)
end

@doc raw"""
	nse_budget_planning!(EP::Model, inputs::Dict, setup::Dict)

Benders planning-problem side of the expected-NSE cap: budgets `vNSEbudget[w] >= 0`, one per
representative period, and `cNSEBudget_planning`: $\sum_w$ `vNSEbudget[w]` $\leq \bar{E}$.
`vNSEbudget[w]` is in the weighted convention (it already carries `Sub_Weights[w]`, i.e.
$\pi_w$), so the cap is an unweighted sum. See [`nse_budget!`](@ref).
"""
function nse_budget_planning!(EP::Model, inputs::Dict, setup::Dict)
    println("NSE Budget Policy Planning Module")

    target = nse_budget_target(inputs, setup)

    if get(setup, :ExpectFeasibleSubproblems, false) == true
        error("NSEBudget = 1 requires ExpectFeasibleSubproblems: false in benders_settings.yml: " *
              "a subproblem whose NSE budget cannot be met is infeasible and must return a feasibility cut.")
    end

    # Weighted convention: vNSEbudget[w] bounds sum_{t in w} omega[t] * NSE[t], and omega[t]
    # carries Sub_Weights[w] (= 8760 * pi_w for scenario-years). Do NOT weight the sum below.
    @variable(EP, vNSEbudget[w = 1:inputs["REP_PERIOD"]]>=0)

    @constraint(EP, cNSEBudget_planning,
        sum(vNSEbudget[w] for w in 1:inputs["REP_PERIOD"])<=target)
end

@doc raw"""
	nse_budget_subperiod!(EP::Model, inputs::Dict, setup::Dict)

Benders subproblem side of the expected-NSE cap, for the representative period
`w = inputs["SubPeriod"]`: `cNSEBudget`: $\sum_{t \in \mathcal{T}_w} \omega_t \sum_{s,z}
\Lambda_{s,t,z} \leq$ `vNSEbudget[w]`. Both sides are in the weighted convention (`inputs["omega"]`
here is the slice of $\omega$ belonging to period $w$, as built by `separate_inputs_subperiods`).
`vNSEbudget[w]` is declared free: it is a linking variable, fixed by `MacroEnergySolvers` at the
value chosen by the planning problem, where its bound lives. See [`nse_budget!`](@ref).
"""
function nse_budget_subperiod!(EP::Model, inputs::Dict, setup::Dict)
    println("NSE Budget Policy Operation Module")

    w = inputs["SubPeriod"]

    @variable(EP, vNSEbudget[[w]])

    # Weighted convention on both sides: omega[t] (= pi_w for scenario-years) multiplies NSE,
    # and vNSEbudget[w] is pi_w * (annual shed budget of w). The subproblem objective carries the
    # same omega[t], so the dual of this row is a price per MWh (times ObjScale).
    @expression(EP, eNSEBudgetExpectedNSE, nse_budget_expected_nse(EP, inputs))

    @constraint(EP, cNSEBudget, eNSEBudgetExpectedNSE<=vNSEbudget[w])
end

# sum_t omega[t] * sum_{s,z} vNSE[s,t,z] over the time steps of `inputs` (the whole horizon in
# the monolithic model, one representative period in a Benders subproblem). Weighted convention.
function nse_budget_expected_nse(EP::Model, inputs::Dict)
    T = inputs["T"]
    Z = inputs["Z"]
    SEG = inputs["SEG"]
    @assert size(EP[:vNSE]) == (SEG, T, Z) "NSE budget: vNSE does not span (SEG, T, Z)"
    # System-wide cap: every zone counts. Zonal prices would need zonal caps.
    return sum(inputs["omega"][t] * EP[:vNSE][s, t, z] for s in 1:SEG, t in 1:T, z in 1:Z)
end

# The cap in model units: NSEBudgetTargetMWh is MWh per year; the model's energy unit is GWh
# when ParameterScale = 1.
function nse_budget_target(inputs::Dict, setup::Dict)
    if !haskey(setup, "NSEBudgetTargetMWh")
        error("NSEBudget is 1 but NSEBudgetTargetMWh (expected annual unserved energy cap, MWh) is not set.")
    end
    target_mwh = setup["NSEBudgetTargetMWh"]
    if !(target_mwh isa Real) || target_mwh < 0
        error("NSEBudgetTargetMWh must be a non-negative number (got $target_mwh).")
    end
    if haskey(inputs, "Weights") && !isapprox(sum(inputs["Weights"]), 8760; rtol = 1e-6) &&
       !haskey(inputs, "SubPeriod")
        @warn "NSE budget: Sub_Weights sum to $(sum(inputs["Weights"])), not 8760; " *
              "NSEBudgetTargetMWh caps sum_t omega[t]*NSE[t], which is then not expected annual energy."
    end
    if nse_budget_removes_voll(setup) && inputs["SEG"] > 1 && !haskey(inputs, "SubPeriod")
        @warn "NSE budget with NSEBudgetRemoveVoLL = 1 and $(inputs["SEG"]) curtailment segments: " *
              "all segments become free and share one cap, so the segment prices no longer order them."
    end
    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1
    return target_mwh / scale_factor
end

# True when the NSE cost term is to be dropped from the objective (budget form).
function nse_budget_removes_voll(setup::Dict)
    return get(setup, "NSEBudget", 0) == 1 && get(setup, "NSEBudgetRemoveVoLL", 1) == 1
end
