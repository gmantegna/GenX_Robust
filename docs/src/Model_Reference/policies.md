# Emission mitigation policies
## Capacity Reserve Margin
```@docs
GenX.cap_reserve_margin!
```

## CO$_2$ Constraint Policy
```@docs
GenX.co2_cap!
```

## Expected Non-Served Energy Budget
```@docs
GenX.nse_budget!
GenX.nse_budget_planning!
GenX.nse_budget_subperiod!
GenX.mes_require_feasible_point!
```

## Energy Share Requirement
```@docs
GenX.load_energy_share_requirement!
GenX.energy_share_requirement!
```

## Minimum Capacity Requirement
```@docs
GenX.minimum_capacity_requirement!
```

## Maximum Capacity Requirement
```@autodocs
Modules = [GenX]
Pages = ["maximum_capacity_requirement.jl"]
```

## Hydrogen Production Demand Requirement (Electrolyzer)
```@docs
GenX.hydrogen_demand!
```

## Hourly clean supply matching constraint
```@docs
GenX.load_hourly_matching_requirement!
GenX.hourly_matching!
```
