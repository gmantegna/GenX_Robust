@doc raw"""
    thermal_plant_effective_capacity(EP::AbstractModel,
                                     inputs::Dict,
                                     resources::Vector{Int},
                                     capres_zone::Int,
                                     timesteps::Vector{Int})::Matrix{Float64})

Effective capacity in a capacity reserve margin zone for certain resources in the given timesteps.
"""
function thermal_plant_effective_capacity(EP,
        inputs,
        resources::Vector{Int},
        capres_zone::Int,
        timesteps::Vector{Int};
        hourly_availability::Bool = false)::Matrix{Float64}
    eff_cap = thermal_plant_effective_capacity.(Ref(EP),
        Ref(inputs),
        resources,
        Ref(capres_zone),
        Ref(timesteps);
        hourly_availability = hourly_availability)
    return reduce(hcat, eff_cap)
end

function thermal_plant_effective_capacity(EP::AbstractModel, inputs::Dict, y, capres_zone::Int;
        hourly_availability::Bool = false)
    T = inputs["T"]
    timesteps = collect(1:T)
    return thermal_plant_effective_capacity(EP, inputs, y, capres_zone, timesteps;
        hourly_availability = hourly_availability)
end

function thermal_plant_effective_capacity(EP::AbstractModel,
        inputs::Dict,
        r_id::Int,
        capres_zone::Int,
        timesteps::Vector{Int};
        hourly_availability::Bool = false)::Vector{Float64}
    y = r_id
    gen = inputs["RESOURCES"]
    capresfactor = derating_factor(gen[y], tag = capres_zone)
    eTotalCap = value.(EP[:eTotalCap][y])

    if hourly_availability # setting CapResThermalHourlyAvailability = 1 (thermal.jl)
        effective_capacity = capresfactor .* inputs["pP_Max"][y, timesteps] .* eTotalCap
    else
        effective_capacity = fill(capresfactor * eTotalCap, length(timesteps))
    end

    if has_maintenance(inputs) && y in ids_with_maintenance(gen)
        adjustment = thermal_maintenance_capacity_reserve_margin_adjustment(EP,
            inputs,
            y,
            capres_zone,
            timesteps)
        effective_capacity = effective_capacity .+ value.(adjustment)
    end

    if y in ids_with(gen, :fusion)
        resource_component = resource_name(gen[y])
        adjustment = thermal_fusion_capacity_reserve_margin_adjustment(
            EP, inputs, resource_component, y, capres_zone, timesteps)
        effective_capacity = effective_capacity .+ value.(adjustment)
    end

    return effective_capacity
end
