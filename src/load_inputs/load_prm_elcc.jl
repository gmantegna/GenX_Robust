@doc raw"""
    load_prm_elcc!(setup::Dict, policies_path::AbstractString, resource_policies_path::AbstractString, inputs::Dict)

Read the inputs of the PRM-with-ELCC-surfaces policy (`PRM_ELCC: 1`; see [`prm_elcc!`](@ref)).

Three files:

- `policies/PRM_ELCC.csv`: header `TRN_MW`, one row: the fixed total reliability need in MW.
- `policies/PRM_ELCC_facets.csv`, long format `Surface,Facet,Axis,Coefficient`: one row per
  (surface, facet, axis). `Axis = intercept` is the facet's constant (MW of ELCC; 0 if absent),
  every other axis row is the facet's slope along that axis (MW of ELCC per MW of axis
  capacity). Facets are named within their surface; a (surface, facet, axis) triple may appear
  once.
- `resources/policy_assignments/Resource_prm_elcc.csv`:
  `Resource,PRM_ELCC_NQC,PRM_ELCC_Axis,PRM_ELCC_Axis_Multiplier`. A firm resource carries an
  NQC (fraction of its power capacity that counts) and a blank axis; a surface resource carries
  an axis name and a multiplier $m_r$ (blank = 1) and no NQC. A resource of the model that is
  missing from the file contributes nothing to the constraint (a warning lists them); a row
  whose resource is not in the model is ignored with a warning. This file is NOT one of GenX's
  generic tagged policy files (`add_policies_to_resources!` ignores it): it is read here and
  matched by resource name.

With `ParameterScale: 1`, `TRN_MW` and the intercepts are divided by `ModelScalingFactor`
(MW to GW); NQCs, slopes and multipliers are dimensionless and untouched.

Stored in `inputs`: `PRM_ELCC_TRN` (model units), `PRM_ELCC_NQC`, `PRM_ELCC_AXIS` (`""` = none),
`PRM_ELCC_MULT` (vectors over resources), `PRM_ELCC_SURFACES` (ordered names),
`PRM_ELCC_FACETS` (surface => ordered facet names), `PRM_ELCC_INTERCEPT` ((surface, facet) =>
constant, model units), `PRM_ELCC_SLOPE` ((surface, facet) => Dict(axis => slope)) and
`PRM_ELCC_AXES` (ordered axis names used by any facet).
"""
function load_prm_elcc!(setup::Dict, policies_path::AbstractString,
        resource_policies_path::AbstractString, inputs::Dict)
    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1.0

    # --- 1. TRN -------------------------------------------------------------------------
    trn_file = "PRM_ELCC.csv"
    trn_df = load_dataframe(joinpath(policies_path, trn_file))
    "TRN_MW" in names(trn_df) || error("$trn_file must have a column TRN_MW (got $(names(trn_df)))")
    nrow(trn_df) == 1 || error("$trn_file must have exactly one row (got $(nrow(trn_df)))")
    trn_mw = trn_df.TRN_MW[1]
    (trn_mw isa Real && !ismissing(trn_mw) && trn_mw >= 0) ||
        error("$trn_file: TRN_MW must be a non-negative number (got $trn_mw)")
    inputs["PRM_ELCC_TRN"] = Float64(trn_mw) / scale_factor

    # --- 2. facets ----------------------------------------------------------------------
    facets_file = "PRM_ELCC_facets.csv"
    fdf = load_dataframe(joinpath(policies_path, facets_file))
    for c in ("Surface", "Facet", "Axis", "Coefficient")
        c in names(fdf) || error("$facets_file must have columns Surface,Facet,Axis,Coefficient (missing $c)")
    end
    surfaces = String[]
    facets = Dict{String, Vector{String}}()
    intercept = Dict{Tuple{String, String}, Float64}()
    slope = Dict{Tuple{String, String}, Dict{String, Float64}}()
    axes = String[]
    seen_rows = Set{Tuple{String, String, String}}()
    for row in eachrow(fdf)
        any(ismissing, (row.Surface, row.Facet, row.Axis, row.Coefficient)) &&
            error("$facets_file: blank cell in row $(NamedTuple(row))")
        s, f, a = strip(string(row.Surface)), strip(string(row.Facet)), strip(string(row.Axis))
        coef = row.Coefficient isa Real ? Float64(row.Coefficient) :
               parse(Float64, strip(string(row.Coefficient)))
        s in surfaces || push!(surfaces, s)
        fl = get!(facets, s, String[])
        if !(f in fl)
            push!(fl, f)
            intercept[(s, f)] = 0.0
            slope[(s, f)] = Dict{String, Float64}()
        end
        key = (s, f, lowercase(a) == "intercept" ? "intercept" : a)
        key in seen_rows && error("$facets_file: duplicate row for surface $s facet $f axis $a")
        push!(seen_rows, key)
        if lowercase(a) == "intercept"
            intercept[(s, f)] = coef / scale_factor
        else
            slope[(s, f)][a] = coef
            a in axes || push!(axes, a)
        end
    end
    isempty(surfaces) &&
        @warn "$facets_file has no facet rows: PRM_ELCC reduces to the firm NQC sum >= TRN."

    # --- 3. resources -------------------------------------------------------------------
    rfile = "Resource_prm_elcc.csv"
    rdf = load_dataframe(joinpath(resource_policies_path, rfile))
    for c in ("Resource", "PRM_ELCC_NQC", "PRM_ELCC_Axis", "PRM_ELCC_Axis_Multiplier")
        c in names(rdf) || error("$rfile must have columns Resource,PRM_ELCC_NQC,PRM_ELCC_Axis," *
                                 "PRM_ELCC_Axis_Multiplier (missing $c)")
    end
    gen = inputs["RESOURCES"]
    G = length(gen)
    rnames = resource_name.(gen)
    nqc = zeros(Float64, G)
    axis = fill("", G)
    mult = zeros(Float64, G)
    seen = falses(G)
    num(x) = ismissing(x) ? nothing :
             x isa Real ? Float64(x) :
             isempty(strip(string(x))) ? nothing : parse(Float64, strip(string(x)))
    str(x) = ismissing(x) ? "" : strip(string(x))
    for row in eachrow(rdf)
        rname = str(row.Resource)
        y = findfirst(==(rname), rnames)
        if y === nothing
            @warn "$rfile: resource $rname is not in the model; row ignored."
            continue
        end
        seen[y] && error("$rfile: resource $rname appears twice")
        seen[y] = true
        q = num(row.PRM_ELCC_NQC)
        a = str(row.PRM_ELCC_Axis)
        m = num(row.PRM_ELCC_Axis_Multiplier)
        if !isempty(a)
            (q === nothing || q == 0) ||
                error("$rfile: resource $rname has both an NQC ($q) and a surface axis ($a); it must be one or the other")
            a in axes || error("$rfile: resource $rname is on axis $a, which no facet of $facets_file uses (axes: $(join(axes, ", ")))")
            axis[y] = a
            mult[y] = m === nothing ? 1.0 : m
        else
            m === nothing || m == 0 || m == 1 ||
                @warn "$rfile: resource $rname has a multiplier ($m) but no axis; the multiplier is ignored."
            nqc[y] = q === nothing ? 0.0 : q
        end
    end
    unseen = rnames[.!seen]
    isempty(unseen) ||
        @warn "PRM_ELCC: $(length(unseen)) resource(s) are not in $rfile and contribute nothing " *
              "to the reliability constraint: $(join(unseen, ", "))"
    unused_axes = [a for a in axes if !(a in axis)]
    isempty(unused_axes) ||
        @warn "PRM_ELCC: axis/axes $(join(unused_axes, ", ")) appear in $facets_file but no resource is assigned to them; their terms are 0."

    inputs["PRM_ELCC_NQC"] = nqc
    inputs["PRM_ELCC_AXIS"] = axis
    inputs["PRM_ELCC_MULT"] = mult
    inputs["PRM_ELCC_SURFACES"] = surfaces
    inputs["PRM_ELCC_FACETS"] = facets
    inputs["PRM_ELCC_INTERCEPT"] = intercept
    inputs["PRM_ELCC_SLOPE"] = slope
    inputs["PRM_ELCC_AXES"] = axes
    println("$trn_file, $facets_file and $rfile Successfully Read!")
    return nothing
end
