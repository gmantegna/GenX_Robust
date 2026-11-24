@doc raw"""
    load_network_data!(setup::Dict, path::AbstractString, inputs_nw::Dict)

Function for reading input parameters related to the electricity transmission network
"""
function load_network_data!(setup::Dict, path::AbstractString, inputs_nw::Dict)
    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    filename = "Network.csv"
    network_var = load_dataframe(joinpath(path, filename))

    as_vector(col::Symbol) = collect(skipmissing(network_var[!, col]))
    to_floats(col::Symbol) = convert(Array{Float64}, as_vector(col))

    # Number of zones in the network
    Z = length(as_vector(:Network_zones))
    inputs_nw["Z"] = Z
    # Number of lines in the network
    L = length(as_vector(:Network_Lines))
    inputs_nw["L"] = L

    # Topology of the network source-sink matrix
    inputs_nw["pNet_Map"] = load_network_map(network_var, Z, L)

    # Transmission capacity of the network (in MW)
    inputs_nw["pTrans_Max"] = to_floats(:Line_Max_Flow_MW) / scale_factor  # convert to GW

    if setup["Trans_Loss_Segments"] == 1
        # Line percentage Loss - valid for case when modeling losses as a fixed percent of absolute value of power flows
        inputs_nw["pPercent_Loss"] = to_floats(:Line_Loss_Percentage)
    elseif setup["Trans_Loss_Segments"] >= 2
        # Transmission line voltage (in kV)
        inputs_nw["kV"] = to_floats(:Line_Voltage_kV)
        # Transmission line resistance (in Ohms) - Used when modeling quadratic transmission losses
        inputs_nw["Ohms"] = to_floats(:Line_Resistance_Ohms)
    end

    ## Inputs for the DC-OPF 
    if setup["DC_OPF"] == 1
        if setup["NetworkExpansion"] == 1
            @warn("Because the DC_OPF flag is active, GenX will not allow any transmission capacity expansion. Set the DC_OPF flag to 0 if you want to optimize tranmission capacity expansion.")
            setup["NetworkExpansion"] = 0
        end
        println("Reading DC-OPF values...")
        # Transmission line voltage (in kV)
        line_voltage_kV = to_floats(:Line_Voltage_kV)
        # Transmission line reactance (in Ohms)
        line_reactance_Ohms = to_floats(:Line_Reactance_Ohms)
        # Line angle limit (in radians)
        inputs_nw["Line_Angle_Limit"] = to_floats(:Angle_Limit_Rad)
        # DC-OPF coefficient for each line (in MW when not scaled, in GW when scaled) 
        # MW = (kV)^2/Ohms 
        inputs_nw["pDC_OPF_coeff"] = ((line_voltage_kV .^ 2) ./ line_reactance_Ohms) /
                                     scale_factor
    end

    # Maximum possible flow after reinforcement for use in linear segments of piecewise approximation
    inputs_nw["pTrans_Max_Possible"] = inputs_nw["pTrans_Max"]

    if setup["NetworkExpansion"] == 1
        # Read between zone network reinforcement costs per peak MW of capacity added
        inputs_nw["pC_Line_Reinforcement"] = to_floats(:Line_Reinforcement_Cost_per_MWyr) /
                                             scale_factor # convert to million $/GW/yr with objective function in millions
        # Maximum reinforcement allowed in MW
        #NOTE: values <0 indicate no expansion possible
        inputs_nw["pMax_Line_Reinforcement"] = map(x -> max(0, x),
            to_floats(:Line_Max_Reinforcement_MW)) / scale_factor # convert to GW
        inputs_nw["pTrans_Max_Possible"] += inputs_nw["pMax_Line_Reinforcement"]
    end

    # Multi-Stage
    if setup["MultiStage"] == 1 
        # Weighted Average Cost of Capital for Transmission Expansion
        if setup["NetworkExpansion"] >= 1
            inputs_nw["transmission_WACC"] = to_floats(:WACC)
            inputs_nw["Capital_Recovery_Period_Trans"] = to_floats(:Capital_Recovery_Period)
        end

        # Max Flow Possible on Each Line
        inputs_nw["pLine_Max_Flow_Possible_MW"] = to_floats(:Line_Max_Flow_Possible_MW) /
                                                  scale_factor # Convert to GW
    end

    # Transmission line (between zone) loss coefficient (resistance/voltage^2)
    inputs_nw["pTrans_Loss_Coef"] = zeros(Float64, L)
    if setup["Trans_Loss_Segments"] == 1
        inputs_nw["pTrans_Loss_Coef"] = inputs_nw["pPercent_Loss"]
    elseif setup["Trans_Loss_Segments"] >= 2
        # If zones are connected, loss coefficient is R/V^2 where R is resistance in Ohms and V is voltage in Volts
        inputs_nw["pTrans_Loss_Coef"] = (inputs_nw["Ohms"] / 10^6) ./
                                        (inputs_nw["kV"] / 10^3)^2 * scale_factor # 1/GW ***
    end

    ## Sets and indices for transmission losses and expansion
    inputs_nw["TRANS_LOSS_SEGS"] = setup["Trans_Loss_Segments"] # Number of segments used in piecewise linear approximations quadratic loss functions
    inputs_nw["LOSS_LINES"] = findall(inputs_nw["pTrans_Loss_Coef"] .!= 0) # Lines for which loss coefficients apply (are non-zero);

    if setup["NetworkExpansion"] == 1
        # Network lines and zones that are expandable have non-negative maximum reinforcement inputs
        inputs_nw["EXPANSION_LINES"] = findall(inputs_nw["pMax_Line_Reinforcement"] .>= 0)
        inputs_nw["NO_EXPANSION_LINES"] = findall(inputs_nw["pMax_Line_Reinforcement"] .< 0)
    end

    load_directional_ratings!(network_var, inputs_nw, setup)

    println(filename * " Successfully Read!")

    load_simultaneous_flow_groups!(path, inputs_nw, setup, scale_factor)

    return network_var
end

@doc raw"""
    load_directional_ratings!(network_var::DataFrame, inputs_nw::Dict, setup::Dict)

Optional `Profile_Forward` / `Profile_Reverse` columns of Network.csv (fractions of the line's
available capacity allowed in the start->end and end->start direction; ported from the chair's
fork). Absent: nothing is stored and no row is built (the model is unchanged). Present: both
columns must be given, one value in [0, 1] per line; a line with a value below 1 must be a loss
line (the directional rows act on `vTAUX_POS`/`vTAUX_NEG`, which exist only for loss lines), and
`Trans_Loss_Segments` must be 1.
"""
function load_directional_ratings!(network_var::DataFrame, inputs_nw::Dict, setup::Dict)
    has_f = hasproperty(network_var, :Profile_Forward)
    has_r = hasproperty(network_var, :Profile_Reverse)
    (has_f || has_r) || return nothing
    (has_f && has_r) ||
        error("Network.csv: Profile_Forward and Profile_Reverse must be given together")
    L = inputs_nw["L"]
    pf = convert(Array{Float64}, collect(skipmissing(network_var[!, :Profile_Forward])))
    pr = convert(Array{Float64}, collect(skipmissing(network_var[!, :Profile_Reverse])))
    (length(pf) == L && length(pr) == L) ||
        error("Network.csv: Profile_Forward/Profile_Reverse need one value per line ($L lines), got $(length(pf))/$(length(pr))")
    for (nm, p) in (("Profile_Forward", pf), ("Profile_Reverse", pr)), l in 1:L
        (isfinite(p[l]) && 0 <= p[l] <= 1) ||
            error("Network.csv: $nm of line $l must be in [0, 1], got $(p[l])")
    end
    setup["Trans_Loss_Segments"] == 1 ||
        error("Network.csv: Profile_Forward/Profile_Reverse are implemented for Trans_Loss_Segments = 1 only")
    for l in 1:L
        if (pf[l] < 1 || pr[l] < 1) && !(l in inputs_nw["LOSS_LINES"])
            error("Network.csv: line $l has a directional rating below 1 but is not a loss line (Line_Loss_Percentage 0), so its directional flow variables do not exist")
        end
    end
    inputs_nw["Profile_Forward"] = pf
    inputs_nw["Profile_Reverse"] = pr
    return nothing
end

@doc raw"""
    load_simultaneous_flow_groups!(path, inputs_nw, setup, scale_factor)

Optional `Simultaneous_Flow_Constraints.csv` in the system folder (ported from the chair's fork).
Columns: `Simultaneous Flow Group`, `Transmission Path` (label only), `Line_Number`, `Direction`
(`forward` = start->end flow `vTAUX_POS`, `reverse` = `vTAUX_NEG`), `limit_MW`. For each group g
and hour t: sum over the group's rows of the directional flow <= limit_g. Absent: nothing stored.
Refused: a direction other than forward/reverse; a `Line_Number` outside 1..L; a line that is not a
loss line; a group with more than one distinct `limit_MW` (the fork takes the first silently); a
(line, direction) listed twice in one group; a negative or non-finite limit;
`Trans_Loss_Segments` other than 1. `limit_MW` is divided by the ParameterScale factor, as
`Line_Max_Flow_MW` is. Stored as `inputs_nw["SimFlowGroups"]`, a vector of named tuples
`(name, lines, forward, limit)` in first-appearance order of the groups.
"""
function load_simultaneous_flow_groups!(path::AbstractString, inputs_nw::Dict, setup::Dict,
        scale_factor)
    filename = "Simultaneous_Flow_Constraints.csv"
    isfile(joinpath(path, filename)) || return nothing
    df = load_dataframe(joinpath(path, filename))
    need = ["Simultaneous Flow Group", "Line_Number", "Direction", "limit_MW"]
    for c in need
        c in names(df) || error("$filename: missing column '$c' (has $(names(df)))")
    end
    setup["Trans_Loss_Segments"] == 1 ||
        error("$filename: simultaneous-flow groups are implemented for Trans_Loss_Segments = 1 only")
    L = inputs_nw["L"]
    LOSS_LINES = inputs_nw["LOSS_LINES"]
    groups = []
    for g in unique(string.(df[!, "Simultaneous Flow Group"]))
        rows = df[string.(df[!, "Simultaneous Flow Group"]) .== g, :]
        limits = unique(Float64.(rows[!, "limit_MW"]))
        length(limits) == 1 ||
            error("$filename: group $g has more than one limit_MW: $limits")
        limit = limits[1]
        (isfinite(limit) && limit >= 0) || error("$filename: group $g has limit_MW $limit")
        lines = Int[]
        fwd = Bool[]
        for row in eachrow(rows)
            d = string(row["Direction"])
            d in ("forward", "reverse") ||
                error("$filename: group $g: direction must be forward or reverse, got $d")
            l = Int(row["Line_Number"])
            1 <= l <= L || error("$filename: group $g: Line_Number $l is outside 1..$L")
            l in LOSS_LINES ||
                error("$filename: group $g: line $l is not a loss line, so its directional flow variables do not exist")
            any(i -> lines[i] == l && fwd[i] == (d == "forward"), eachindex(lines)) &&
                error("$filename: group $g lists line $l $d twice")
            push!(lines, l)
            push!(fwd, d == "forward")
        end
        push!(groups, (name = g, lines = lines, forward = fwd, limit = limit / scale_factor))
    end
    inputs_nw["SimFlowGroups"] = groups
    println(filename * " Successfully Read! ($(length(groups)) groups)")
    return nothing
end

@doc raw"""
    load_network_map_from_list(network_var::DataFrame, Z, L, list_columns)

Loads the network map from a list-style interface
```
..., Network_Lines, Start_Zone, End_Zone, ...
                 1,           1,       2,
                 2,           1,       3,
```
"""
function load_network_map_from_list(network_var::DataFrame, Z, L, list_columns)
    start_col, end_col = list_columns
    mat = zeros(L, Z)
    start_zones = collect(skipmissing(network_var[!, start_col]))
    end_zones = collect(skipmissing(network_var[!, end_col]))
    for l in 1:L
        mat[l, start_zones[l]] = 1
        mat[l, end_zones[l]] = -1
    end
    mat
end

@doc raw"""
    load_network_map_from_matrix(network_var::DataFrame, Z, L)

Loads the network map from a matrix-style interface
```
..., Network_Lines, z1, z2, z3, ...
                 1,  1, -1,  0,
                 2,  1,  0, -1,
```
This is equivalent to the list-style interface where the zone zN with entry +1 is the
starting zone of the line and the zone with entry -1 is the ending zone of the line.
"""
function load_network_map_from_matrix(network_var::DataFrame, Z, L)
    # Topology of the network source-sink matrix
    network_map_matrix_format_deprecation_warning()
    col = findall(s -> s == "z1", names(network_var))[1]
    mat = Matrix{Float64}(network_var[1:L, col:(col + Z - 1)])
end

function load_network_map(network_var::DataFrame, Z, L)
    columns = names(network_var)

    list_columns = ["Start_Zone", "End_Zone"]
    has_network_list = all([c in columns for c in list_columns])

    zones_as_strings = ["z" * string(i) for i in 1:Z]
    has_network_matrix = all([c in columns for c in zones_as_strings])

    instructions = """The transmission network should be specified in the form of a matrix
           (with columns z1, z2, ... zN) or in the form of lists (with Start_Zone, End_Zone),
           but not both. See the documentation for examples."""

    if has_network_list && has_network_matrix
        error("two types of transmission network map were provided.\n" * instructions)
    elseif !(has_network_list || has_network_matrix)
        error("no transmission network map was detected.\n" * instructions)
    elseif has_network_list
        load_network_map_from_list(network_var, Z, L, list_columns)
    elseif has_network_matrix
        load_network_map_from_matrix(network_var, Z, L)
    end
end

function network_map_matrix_format_deprecation_warning()
    @warn """Specifying the network map as a matrix is deprecated as of v0.4
  and will be removed in v0.5. Instead, use the more compact list-style format.

  ..., Network_Lines, Start_Zone, End_Zone, ...
                   1,          1,        2,
                   2,          1,        3,
                   3,          2,        3,
  """ maxlog=1
end
