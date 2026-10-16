@doc raw"""
	prm_elcc!(EP::Model, inputs::Dict, setup::Dict)

Planning reserve margin with ELCC accreditation, in the RESOLVE form (a fixed total reliability
need met by firm capacity at a net qualifying capacity plus piecewise-linear ELCC surfaces).
Activated by `PRM_ELCC: 1` in `genx_settings.yml`; inputs in [`load_prm_elcc!`](@ref). One
annual, system-wide constraint (zone-agnostic: every zone's capacity counts towards the one need).

Let $\mathcal{F}$ be the firm resources (those with an NQC $q_r$), $\mathcal{S}$ the surfaces,
$\mathcal{P}_s$ the facets of surface $s$, $\mathcal{A}$ the axes, $\mathcal{G}_a$ the
resources assigned to axis $a$ with multipliers $m_r$, and $\Delta^{\text{total}}_r$ the total
**power** capacity of resource $r$ (`eTotalCap`; storage counts by its discharge capacity, never
by energy). With a free variable $\mathrm{ELCC}_s$ (`vELCC[s]`) per surface,

```math
\begin{aligned}
	\sum_{r \in \mathcal{F}} q_r \, \Delta^{\text{total}}_r + \sum_{s \in \mathcal{S}} \mathrm{ELCC}_s &\geq \mathrm{TRN}
	& \text{(cPRM\_ELCC)} \\
	\mathrm{ELCC}_s &\leq b_{s,f} + \sum_{a \in \mathcal{A}} c_{s,f,a} \sum_{r \in \mathcal{G}_a} m_r \, \Delta^{\text{total}}_r
	& \forall s \in \mathcal{S},\ f \in \mathcal{P}_s \quad \text{(cELCCFacet}[s,f])
\end{aligned}
```

so that $\mathrm{ELCC}_s$ is at most the lower envelope of the facets, i.e. the concave
piecewise-linear surface $\min_f (b_{s,f} + c_{s,f}^\top x)$, which is exactly RESOLVE's ELCC
surface constraint (a linearization of the ELCC of the portfolio, computed by a separate
reliability model). In an optimum where the reliability need binds, $\mathrm{ELCC}_s$ sits on
the lowest facet at the portfolio and the dual of `cPRM_ELCC` is the capacity price in
\$/MW-yr (times `ModelScalingFactor` when `ParameterScale: 1`, i.e. the raw dual is in M\$/GW).

The intercepts and TRN are in MW (GW under `ParameterScale: 1`); the slopes $c_{s,f,a}$ are MW
of ELCC per MW of axis capacity and the multipliers $m_r$ are dimensionless, so the constraint
needs no further scaling. Resources missing from the assignment file have $q_r = 0$ and no axis.

Expressions written for the output (`prm_elcc.csv`, [`write_prm_elcc`](@ref)):
`ePRM_ELCC_Firm` ($\sum_{\mathcal F} q_r \Delta_r$), `ePRM_ELCC_AxisCap[a]`
($\sum_{\mathcal G_a} m_r \Delta_r$) and `ePRM_ELCC_FacetRHS[s,f]` (the right-hand side of each
facet).

Two things to know when writing the inputs:

- An axis is a global name: a facet row of ANY surface that names axis `a` takes the same
  $\sum_{r \in \mathcal{G}_a} m_r \Delta^{\text{total}}_r$, so a resource's capacity counts in
  every surface that uses its axis. Nothing warns about it: if two surfaces share an axis name the
  same MW are accredited twice (once per surface). Give each resource group its own axis name
  unless that double count is intended.
- $\Delta^{\text{total}}_r$ is `eTotalCap[r]`, the one total-capacity expression of the resource.
  For a co-located `VRE_STOR` resource it is the inverter (grid-connection) MW, not the DC solar,
  wind or storage components, which have their own expressions (`eTotalCap_SOLAR`, ...). An NQC or
  axis multiplier on such a resource applies to its inverter MW.

This is a planning-level constraint on capacities only; it is implemented for the monolithic
model. With `Benders: 1` the model build stops with an error (it would belong in the planning
problem, where it is a constraint on linking variables only, but that path is not implemented).
"""
function prm_elcc!(EP::Model, inputs::Dict, setup::Dict)
    println("PRM with ELCC Surfaces Policy Module")

    G = inputs["G"]
    nqc = inputs["PRM_ELCC_NQC"]
    axis = inputs["PRM_ELCC_AXIS"]
    mult = inputs["PRM_ELCC_MULT"]
    surfaces = inputs["PRM_ELCC_SURFACES"]
    facets = inputs["PRM_ELCC_FACETS"]
    intercept = inputs["PRM_ELCC_INTERCEPT"]
    slope = inputs["PRM_ELCC_SLOPE"]
    axes = inputs["PRM_ELCC_AXES"]
    trn = inputs["PRM_ELCC_TRN"]

    # eTotalCap[y] is an AffExpr for buildable/retirable resources and a plain number for a
    # fixed one; add_to_expression! handles both, so build the sums explicitly.
    firm = AffExpr(0.0)
    for y in 1:G
        nqc[y] == 0 || add_to_expression!(firm, nqc[y], EP[:eTotalCap][y])
    end
    EP[:ePRM_ELCC_Firm] = firm

    axis_cap = Dict{String, AffExpr}()
    for a in axes
        ex = AffExpr(0.0)
        for y in 1:G
            axis[y] == a && add_to_expression!(ex, mult[y], EP[:eTotalCap][y])
        end
        axis_cap[a] = ex
    end
    EP[:ePRM_ELCC_AxisCap] = axis_cap

    # One free ELCC variable per surface: its only bounds are the facets below.
    @variable(EP, vELCC[s in surfaces])

    facet_rhs = Dict{Tuple{String, String}, AffExpr}()
    for s in surfaces, f in facets[s]
        rhs = AffExpr(intercept[(s, f)])
        for (a, c) in slope[(s, f)]
            add_to_expression!(rhs, c, axis_cap[a])
        end
        facet_rhs[(s, f)] = rhs
    end
    EP[:ePRM_ELCC_FacetRHS] = facet_rhs

    # ELCC_s <= every facet of its surface (the lower envelope).
    @constraint(EP, cELCCFacet[s in surfaces, f in facets[s]],
        vELCC[s]<=facet_rhs[(s, f)])

    # The reliability need: firm NQC-accredited capacity plus the surfaces' ELCC.
    accredited = copy(firm)
    for s in surfaces
        add_to_expression!(accredited, vELCC[s])
    end
    EP[:ePRM_ELCC_Accredited] = accredited
    @constraint(EP, cPRM_ELCC, accredited>=trn)

    return nothing
end
