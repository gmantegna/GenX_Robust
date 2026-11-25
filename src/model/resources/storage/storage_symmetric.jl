@doc raw"""
	storage_symmetric!(EP::AbstractModel, inputs::Dict, setup::Dict)

Sets up variables and constraints specific to storage resources with symmetric charge and discharge capacities. See ```storage()``` in ```storage.jl``` for description of constraints.
"""
function storage_symmetric!(EP::AbstractModel, inputs::Dict, setup::Dict)
    # Set up additional variables, constraints, and expressions associated with storage resources with symmetric charge & discharge capacity
    # (e.g. most electrochemical batteries that use same components for charge & discharge)
    # STOR = 1 corresponds to storage with distinct power and energy capacity decisions but symmetric charge/discharge power ratings

    println("Storage Resources with Symmetric Charge/Discharge Capacity Module")

    OperationalReserves = setup["OperationalReserves"]
    CapacityReserveMargin = setup["CapacityReserveMargin"]

    T = inputs["T"]     # Number of time steps (hours)

    # operational constraints apply only to the operational storage set (op-group
    # representatives + standalone), sized by the group capacity expression
    STOR_SYMMETRIC = intersect(inputs["STOR_SYMMETRIC"], inputs["STOR_OPERATIONAL"])

    ### Constraints ###

    # Storage discharge and charge power (and reserve contribution) related constraints for symmetric storage resources:
    if OperationalReserves == 1
        storage_symmetric_operational_reserves!(EP, inputs, setup)
    else
        if CapacityReserveMargin > 0
            @constraints(EP,
                begin
                    # Maximum charging rate (including virtual charging to move energy held in reserve back to available storage) must be less than symmetric power rating
                    # Max simultaneous charge and discharge cannot be greater than capacity
                    [y in STOR_SYMMETRIC, t in 1:T],
                    EP[:vP][y, t] + EP[:vCHARGE][y, t] + EP[:vCAPRES_discharge][y, t] +
                    EP[:vCAPRES_charge][y, t] <= EP[:eTotalCapStorOp][y]
                end)
        else
            @constraints(EP,
                begin
                    # Maximum charging rate (including virtual charging to move energy held in reserve back to available storage) must be less than symmetric power rating
                    # Max simultaneous charge and discharge cannot be greater than capacity
                    [y in STOR_SYMMETRIC, t in 1:T],
                    EP[:vP][y, t] + EP[:vCHARGE][y, t] <= EP[:eTotalCapStorOp][y]
                end)
        end
    end

    # Simultaneous charge and discharge within the hour, off unless StorSimultaneousCharge == 1:
    # vP + vCHARGE <= 0.5 x (pP_Max + 1) x capacity, on the resources whose discharge profile pP_Max
    # is below 1 in some hour (charging profile 1); elsewhere the row equals the one above
    if get(setup, "StorSimultaneousCharge", 0) == 1
        STOR_SIM = [y for y in STOR_SYMMETRIC if any(inputs["pP_Max"][y, t] < 1.0 for t in 1:T)]
        @constraint(EP, cStorSimCharge[y in STOR_SIM, t in 1:T],
            EP[:vP][y, t] + EP[:vCHARGE][y, t] <=
            0.5 * (inputs["pP_Max"][y, t] + 1.0) * EP[:eTotalCapStorOp][y])
    end
end

@doc raw"""
	storage_symmetric_operational_reserves!(EP::AbstractModel, inputs::Dict)

Sets up variables and constraints specific to storage resources with symmetric charge and discharge capacities when reserves are modeled. See ```storage()``` in ```storage.jl``` for description of constraints.
"""
function storage_symmetric_operational_reserves!(EP::AbstractModel, inputs::Dict, setup::Dict)
    T = inputs["T"]
    CapacityReserveMargin = setup["CapacityReserveMargin"] > 0

    SYMMETRIC = intersect(inputs["STOR_SYMMETRIC"], inputs["STOR_OPERATIONAL"])

    REG = intersect(SYMMETRIC, inputs["REG"])
    RSV = intersect(SYMMETRIC, inputs["RSV"])

    vP = EP[:vP]
    vCHARGE = EP[:vCHARGE]
    vREG_charge = EP[:vREG_charge]
    vRSV_charge = EP[:vRSV_charge]
    vREG_discharge = EP[:vREG_discharge]
    vRSV_discharge = EP[:vRSV_discharge]
    eTotalCap = EP[:eTotalCapStorOp]

    # Maximum charging rate plus contribution to regulation down must be less than symmetric power rating
    # Max simultaneous charge and discharge rates cannot be greater than symmetric charge/discharge capacity
    expr = @expression(EP, [y in SYMMETRIC, t in 1:T], vP[y, t]+vCHARGE[y, t])
    add_similar_to_expression!(expr[REG, :], vREG_charge[REG, :])
    add_similar_to_expression!(expr[REG, :], vREG_discharge[REG, :])
    add_similar_to_expression!(expr[RSV, :], vRSV_discharge[RSV, :])
    if CapacityReserveMargin
        vCAPRES_charge = EP[:vCAPRES_charge]
        vCAPRES_discharge = EP[:vCAPRES_discharge]
        add_similar_to_expression!(expr[SYMMETRIC, :], vCAPRES_charge[SYMMETRIC, :])
        add_similar_to_expression!(expr[SYMMETRIC, :], vCAPRES_discharge[SYMMETRIC, :])
    end
    @constraint(EP, [y in SYMMETRIC, t in 1:T], expr[y, t]<=eTotalCap[y])
end
