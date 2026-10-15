@doc raw"""
	write_prm_elcc(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)

Write `prm_elcc.csv` for the PRM-with-ELCC-surfaces policy ([`prm_elcc!`](@ref)): long format,
columns `Item,Surface,Facet,Axis,Resource,Value,Unit`, one row per item:

| Item | keyed by | Value |
|---|---|---|
| `TRN` | | the reliability need, MW |
| `Dual` | | the dual of `cPRM_ELCC`: the capacity price, \$/MW-yr (0 when the need does not bind) |
| `TotalAccredited` | | firm accreditation + sum of the surfaces' ELCC, MW |
| `FirmAccredited` | Resource | $q_r \Delta_r$ for each firm resource, MW |
| `AxisCapacity` | Axis | $\sum_{r \in \mathcal{G}_a} m_r \Delta_r$, MW |
| `SurfaceELCC` | Surface | `vELCC[s]`, MW |
| `FacetValue` | Surface, Facet | the facet's right-hand side at the solution, MW |
| `FacetSlack` | Surface, Facet | right-hand side minus `vELCC[s]`, MW (0 = tight) |
| `FacetDual` | Surface, Facet | the facet's dual, \$/MW-yr |
| `BindingFacet` | Surface, Facet | the facet with the smallest slack (the lowest plane at the portfolio); Value = its slack |

Capacities and MW quantities are multiplied by `ModelScalingFactor` when `ParameterScale: 1`
(GW to MW); duals are multiplied by it too (M\$/GW to \$/MW-yr), as `write_minimum_capacity_requirement`
does. Dual rows are written only when the model has duals.
"""
function write_prm_elcc(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1.0
    gen = inputs["RESOURCES"]
    G = inputs["G"]
    nqc = inputs["PRM_ELCC_NQC"]
    surfaces = inputs["PRM_ELCC_SURFACES"]
    facets = inputs["PRM_ELCC_FACETS"]
    axes = inputs["PRM_ELCC_AXES"]
    duals = has_duals(EP)

    df = DataFrame(Item = String[], Surface = Union{Missing, String}[],
        Facet = Union{Missing, String}[], Axis = Union{Missing, String}[],
        Resource = Union{Missing, String}[], Value = Float64[], Unit = String[])
    add!(item, val, unit; surface = missing, facet = missing, axis = missing, resource = missing) =
        push!(df, (item, surface, facet, axis, resource, val, unit))

    add!("TRN", inputs["PRM_ELCC_TRN"] * scale_factor, "MW")
    if duals
        add!("Dual", dual(EP[:cPRM_ELCC]) * scale_factor, "\$/MW-yr")
    end
    firm = value(EP[:ePRM_ELCC_Firm])
    elcc = Dict(s => value(EP[:vELCC][s]) for s in surfaces)
    add!("TotalAccredited", (firm + sum(values(elcc); init = 0.0)) * scale_factor, "MW")
    for y in 1:G
        nqc[y] == 0 && continue
        add!("FirmAccredited", nqc[y] * value(EP[:eTotalCap][y]) * scale_factor, "MW";
            resource = resource_name(gen[y]))
    end
    for a in axes
        add!("AxisCapacity", value(EP[:ePRM_ELCC_AxisCap][a]) * scale_factor, "MW"; axis = a)
    end
    for s in surfaces
        add!("SurfaceELCC", elcc[s] * scale_factor, "MW"; surface = s)
        slacks = Dict{String, Float64}()
        for f in facets[s]
            rhs = value(EP[:ePRM_ELCC_FacetRHS][(s, f)])
            slacks[f] = rhs - elcc[s]
            add!("FacetValue", rhs * scale_factor, "MW"; surface = s, facet = f)
            add!("FacetSlack", slacks[f] * scale_factor, "MW"; surface = s, facet = f)
            if duals
                add!("FacetDual", dual(EP[:cELCCFacet][s, f]) * scale_factor, "\$/MW-yr";
                    surface = s, facet = f)
            end
        end
        if !isempty(facets[s])
            fb = argmin(f -> slacks[f], facets[s])
            add!("BindingFacet", slacks[fb] * scale_factor, "MW"; surface = s, facet = fb)
        end
    end
    CSV.write(joinpath(path, "prm_elcc.csv"), df)
    return nothing
end
