@doc raw"""
    load_co2_cap!(setup::Dict, path::AbstractString, inputs::Dict)

Read input parameters related to CO$_2$ emissions cap constraints
"""
function load_co2_cap!(setup::Dict, path::AbstractString, inputs::Dict)
    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    filename = "CO2_cap_slack.csv"
    if isfile(joinpath(path, filename))
        df = load_dataframe(joinpath(path, filename))
        inputs["dfCO2Cap_slack"] = df
        inputs["dfCO2Cap_slack"][!, :PriceCap] ./= scale_factor # Million $/kton if scaled, $/ton if not scaled
    end

    filename = "CO2_cap.csv"
    df = load_dataframe(joinpath(path, filename))

    inputs["dfCO2Cap"] = df
    mat = extract_matrix_from_dataframe(df, "CO_2_Cap_Zone")
    inputs["dfCO2CapZones"] = mat
    inputs["NCO2Cap"] = size(mat, 2)

    # Emission limits
    if setup["CO2Cap"] == 1
        #  CO2 emissions cap in mass
        # note the default inputs is in million tons
        # when scaled, the constraint unit is kton
        # when not scaled, the constraint unit is ton
        mat = extract_matrix_from_dataframe(df, "CO_2_Max_Mtons")
        inputs["dfMaxCO2"] = mat * 1e6 / scale_factor

    elseif setup["CO2Cap"] == 2 || setup["CO2Cap"] == 3
        #  CO2 emissions rate applied per MWh
        mat = extract_matrix_from_dataframe(df, "CO_2_Max_tons_MWh")
        # no scale_factor is needed since this is a ratio
        inputs["dfMaxCO2Rate"] = mat
    end

    println(filename * " Successfully Read!")
end

@doc raw"""
    load_co2_import_factors!(network_var::DataFrame, inputs::Dict, setup::Dict)

Optional import-CO2 factors of lines (ported from the chair's fork d58adffd,
`load_co2_cap_trans!`): columns `EF_Forward_CO_2_Cap_Zone_c` and `EF_Reverse_CO_2_Cap_Zone_c` of
`Network.csv` (t/MWh of start->end and end->start flow counted against CO2 cap `c`), one column
per cap of `CO2_cap.csv`. Absent: nothing is stored and the CO2 rows are unchanged. Present: both
prefixes, one finite value per line and cap; if every factor is 0 nothing is stored; otherwise stored as `inputs["dfTransEFForward"]`,
`inputs["dfTransEFReverse"]` (L x NCO2Cap). Refused: one prefix without the other; a column count
other than the number of caps; a nonzero factor with `CO2Cap` other than 1 (the term exists for the
mass cap only); a nonzero forward factor on a line that is not a loss line (the term reads
`vTAUX_POS`/`vTAUX_NEG`, which exist only for loss lines); `Trans_Loss_Segments` other than 1;
`ParameterScale: 1`. A line with a zero forward and a nonzero reverse factor for a cap contributes
nothing to that cap (the fork's filter on the forward factor, kept; a warning names the lines).
"""
function load_co2_import_factors!(network_var::DataFrame, inputs::Dict, setup::Dict)
    PF, PR = "EF_Forward_CO_2_Cap_Zone", "EF_Reverse_CO_2_Cap_Zone"
    has_f = !isempty(find_matrix_columns_in_dataframe(network_var, PF))
    has_r = !isempty(find_matrix_columns_in_dataframe(network_var, PR))
    (has_f || has_r) || return nothing
    (has_f && has_r) ||
        error("Network.csv: $(PF)_* and $(PR)_* columns must be given together")
    L, NC = inputs["L"], inputs["NCO2Cap"]
    ef = extract_matrix_from_dataframe(network_var, PF)
    er = extract_matrix_from_dataframe(network_var, PR)
    (size(ef) == (L, NC) && size(er) == (L, NC)) ||
        error("Network.csv: import-CO2 factors need one value per line ($L lines) and one column " *
              "per CO2 cap ($NC caps); got $(size(ef)) forward and $(size(er)) reverse")
    ef = convert(Matrix{Float64}, ef)
    er = convert(Matrix{Float64}, er)
    (all(isfinite, ef) && all(isfinite, er)) ||
        error("Network.csv: import-CO2 factors must be finite")
    any_nonzero = any(!=(0), ef) || any(!=(0), er)
    (any_nonzero && setup["CO2Cap"] != 1) &&
        error("Network.csv: import-CO2 factors are implemented for the mass cap (CO2Cap: 1) only; " *
              "got CO2Cap: $(setup["CO2Cap"]) with a nonzero factor")
    any_nonzero || return nothing    # all factors 0: nothing stored, the CO2 rows are unchanged
    setup["Trans_Loss_Segments"] == 1 ||
        error("Network.csv: import-CO2 factors are implemented for Trans_Loss_Segments = 1 only")
    setup["ParameterScale"] == 1 &&
        error("Network.csv: import-CO2 factors are not supported with ParameterScale: 1")
    for cap in 1:NC, l in 1:L
        if ef[l, cap] != 0 && !(l in inputs["LOSS_LINES"])
            error("Network.csv: line $l has a nonzero $(PF)_$cap but is not a loss line " *
                  "(Line_Loss_Percentage 0), so its directional flow variables do not exist")
        end
    end
    rev_only = [(l, cap) for cap in 1:NC, l in 1:L if ef[l, cap] == 0 && er[l, cap] != 0]
    isempty(rev_only) ||
        @warn "Network.csv: (line, cap) $(rev_only) have a reverse import-CO2 factor but a zero " *
              "forward factor; as in the chair's fork, such a line adds nothing to that cap"
    inputs["dfTransEFForward"] = ef
    inputs["dfTransEFReverse"] = er
    return nothing
end
