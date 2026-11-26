@doc raw"""
	storage_symmetric!(EP::Model, inputs::Dict, setup::Dict)

Sets up constraints specific to storage resources with symmetric charge and discharge capacities. See ```storage()``` in ```storage.jl``` for description of constraints.
"""
function storage_symmetric!(EP::Model, inputs::Dict, setup::Dict)
    # Set up additional constraints associated with storage resources with symmetric charge & discharge capacity
    # (e.g. most electrochemical batteries that use same components for charge & discharge)
    # STOR = 1 corresponds to storage with distinct power and energy capacity decisions but symmetric charge/discharge power ratings

    println("Storage Resources with Symmetric Charge/Discharge Capacity Module")

    OperationalReserves = setup["OperationalReserves"] == 1
    CapacityReserveMargin = setup["CapacityReserveMargin"] > 0

    T = inputs["T"]     # Number of time steps (hours)

    SYMMETRIC = inputs["STOR_SYMMETRIC"]
    eTotalCap = stor_cap_power(EP)   # group capacity if operational groups exist
    vP = EP[:vP]
    vCHARGE = EP[:vCHARGE]

    ### Constraints ###

    # Maximum charging rate plus contribution to regulation down must be less than symmetric power rating
    # Max simultaneous charge and discharge rates cannot be greater than symmetric charge/discharge capacity
    expr = @expression(EP, [y in SYMMETRIC, t in 1:T], vP[y, t] + vCHARGE[y, t])

    if OperationalReserves
        REG = intersect(SYMMETRIC, inputs["REG"])
        RSV = intersect(SYMMETRIC, inputs["RSV"])
        vREG_charge = EP[:vREG_charge]
        vREG_discharge = EP[:vREG_discharge]
        vRSV_discharge = EP[:vRSV_discharge]
        add_similar_to_expression!(expr[REG, :], vREG_charge[REG, :])
        add_similar_to_expression!(expr[REG, :], vREG_discharge[REG, :])
        add_similar_to_expression!(expr[RSV, :], vRSV_discharge[RSV, :])
    end

    if CapacityReserveMargin
        # Maximum charging rate (including virtual charging to move energy held in reserve back to available storage) must be less than symmetric power rating
        vCAPRES_charge = EP[:vCAPRES_charge]
        vCAPRES_discharge = EP[:vCAPRES_discharge]
        add_similar_to_expression!(expr[SYMMETRIC, :], vCAPRES_charge[SYMMETRIC, :])
        add_similar_to_expression!(expr[SYMMETRIC, :], vCAPRES_discharge[SYMMETRIC, :])
    end

    @constraint(EP, [y in SYMMETRIC, t in 1:T], expr[y, t]<=eTotalCap[y])

    # Simultaneous charge and discharge within the hour (StorSimultaneousCharge: 1, ported from the
    # lab's rung-B branch base/d58-rungb f63283ac): for the operational symmetric storage (a group's
    # representative or ungrouped storage) whose pP_Max is below 1 in some hour,
    # vP + vCHARGE <= 0.5 x (pP_Max + 1) x capacity (group capacity for a representative), every hour.
    # Reserve and virtual charge/discharge terms are not in this row (as on f63283ac).
    if get(setup, "StorSimultaneousCharge", 0) == 1
        STOR_SIM = [y for y in intersect(SYMMETRIC, get(inputs, "STOR_OPERATIONAL", SYMMETRIC))
                    if any(inputs["pP_Max"][y, t] < 1.0 for t in 1:T)]
        @constraint(EP, cStorSimCharge[y in STOR_SIM, t in 1:T],
            vP[y, t] + vCHARGE[y, t] <= 0.5 * (inputs["pP_Max"][y, t] + 1.0) * eTotalCap[y])
    end
end
