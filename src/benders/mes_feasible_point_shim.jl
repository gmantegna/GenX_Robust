const MES_FEASIBLE_POINT_SHIM_APPLIED = Ref(false)

@doc raw"""
	mes_require_feasible_point!()

Runtime shim for a defect of `MacroEnergySolvers` 0.2.2 that makes its feasibility-cut path
unreliable with HiGHS. **Opt-in; nothing in GenX calls it.** Call it once per process
(`@everywhere GenX.mes_require_feasible_point!()` with distributed subproblems) before
`MacroEnergySolvers.benders`.

The defect: `MacroEnergySolvers.solve_subproblem` decides "the subproblem solved" by
`has_values(m)`, i.e. `primal_status(m) != NO_SOLUTION`. HiGHS simplex, when it proves
infeasibility from a warm-started basis (every Benders iteration after the first), returns
`termination_status = INFEASIBLE` with `primal_status = INFEASIBLE_POINT`, so `has_values` is
true. MES then takes the objective of that infeasible point and the duals of the fixing rows
(an infeasibility certificate, not a subgradient) and adds them as an *optimality* cut. The cut
is invalid: it cuts off the optimum, and Benders stops with `NEGATIVE GAP` and a lower bound
above the true optimal value. On the first iteration HiGHS presolve usually detects the
infeasibility, returns `NO_SOLUTION`, and the feasibility cut is generated correctly, which is
why the defect is easy to miss. It matters for any model whose subproblems can be infeasible,
e.g. the expected-NSE budget ([`nse_budget!`](@ref)).

The shim redefines that one method with the test `primal_status(m) == FEASIBLE_POINT`
in place of `has_values(m)`; everything else is MES's code. Remove it once MES is fixed.
"""
function mes_require_feasible_point!()
    MES_FEASIBLE_POINT_SHIM_APPLIED[] && return nothing
    @warn "Overriding MacroEnergySolvers.solve_subproblem: a subproblem counts as solved only " *
          "if primal_status == FEASIBLE_POINT (MES 0.2.2 tests has_values)."
    @eval MacroEnergySolvers function solve_subproblem(m::Model, planning_sol::NamedTuple,
            linking_variables_sub::Vector{String}, expect_feasible_subproblems::Bool)
        fix_linking_variables!(m, planning_sol, linking_variables_sub)
        optimize!(m)
        if primal_status(m) == JuMP.MOI.FEASIBLE_POINT
            op_cost = objective_value(m)
            lambda = [dual(FixRef(variable_by_name(m, y))) for y in linking_variables_sub]
            theta_coeff = 1
        elseif expect_feasible_subproblems == true
            error("The subproblem did not return a feasible point (termination status " *
                  "$(termination_status(m))), but ExpectFeasibleSubproblems = true. " *
                  "Set it to false to generate feasibility cuts.")
        else
            @info "Subproblem is infeasible, generating feasibility cut..."
            # Phase-I problem of MES: min slack_max, every affine row relaxed by slack_max
            unfix.(m[:slack_max])
            objfun = objective_function(m)
            @objective(m, Min, m[:slack_max])
            optimize!(m)
            if primal_status(m) != JuMP.MOI.FEASIBLE_POINT
                error("Feasibility subproblem has no feasible point (termination status " *
                      "$(termination_status(m))); this should not happen. Check the model.")
            end
            op_cost = objective_value(m)
            lambda = [dual(FixRef(variable_by_name(m, y))) for y in linking_variables_sub]
            theta_coeff = 0
            fix.(m[:slack_max], 0.0)
            @objective(m, Min, objfun)
        end
        return (op_cost = op_cost, lambda = lambda, theta_coeff = theta_coeff)
    end
    MES_FEASIBLE_POINT_SHIM_APPLIED[] = true
    return nothing
end
