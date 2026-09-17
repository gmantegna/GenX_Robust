# Load-time patch for a defect of MacroEnergySolvers (MES) 0.2.2. See `mes_require_feasible_point!`.

"True once `MacroEnergySolvers.solve_subproblem` has been replaced in this process."
const MES_FEASIBLE_POINT_SHIM_APPLIED = Ref(false)

"The only MES version whose `solve_subproblem` the patch is a copy of."
const MES_PATCHED_VERSION = v"0.2.2"

"Set this environment variable to `1` before `using GenX` to keep stock MES (for testing)."
const MES_PATCH_DISABLE_ENV = "GENX_DISABLE_MES_PATCH"

@doc raw"""
	mes_require_feasible_point!()

Patch for a defect of `MacroEnergySolvers` (MES) 0.2.2 that makes its feasibility-cut path
wrong with HiGHS. **Applied automatically when GenX is loaded** (`GenX.__init__`), on every
process that loads GenX, hence on every Distributed worker. Because the replacement happens at
load time it is visible to every call made after `using GenX`, from the top level or from
inside a function (no world-age constraint). One `@info` line reports it.

- Escape hatch: `ENV["GENX_DISABLE_MES_PATCH"] = "1"` before `using GenX` keeps stock MES.
- Version guard: only MES 0.2.2 is patched. With any other MES version nothing is done; if the
  source of that version still shows the defective branch a warning is issued.
- Not applied while another package that depends on GenX is being precompiled (evaluating into
  MES is not allowed then); it is applied when that package is loaded.
- `GenX.MES_FEASIBLE_POINT_SHIM_APPLIED[]` tells whether the replacement was made.

Calling this function by hand is never needed; it is kept as an idempotent public guard. It
applies the patch if it is not applied yet (escape hatch set, then a change of mind) and
returns whether the patch is in place. Only in that late case does the world-age rule apply:
the replacement is visible to code that starts running after the call returns to top level,
or that is reached through `Base.invokelatest`; GenX's own Benders runner calls MES through
`Base.invokelatest` for that reason.

**The defect.** `MacroEnergySolvers.solve_subproblem` decides "the subproblem solved" by
`has_values(m)`, i.e. `primal_status(m) != NO_SOLUTION`. HiGHS simplex, when it proves
infeasibility from a warm-started basis (every Benders iteration after the first), returns
`termination_status = INFEASIBLE` with `primal_status = INFEASIBLE_POINT`, so `has_values` is
true. MES then takes the objective of that infeasible point and the duals of the fixing rows
(an infeasibility certificate, not a subgradient) and adds them as an *optimality* cut. The cut
is invalid: it cuts off the optimum, the upper bound is corrupted too, and Benders stops with
`NEGATIVE GAP` and a lower bound above the true optimal value. On the first iteration HiGHS
presolve usually detects the infeasibility, returns `NO_SOLUTION`, and the feasibility cut is
generated correctly, which is why the defect is easy to miss. It matters for any model whose
subproblems can be infeasible, e.g. the expected-NSE budget ([`nse_budget!`](@ref)).
`MacroEnergySolvers.solve_planning_problem` has the same `has_values` test; it is NOT patched
(an infeasible master that returns a point has not been observed).

**The replacement** is MES 0.2.2's method with three differences:

1. *Status test.* The subproblem counts as solved, and yields an optimality cut, only if
   `primal_status` is `FEASIBLE_POINT` or `NEARLY_FEASIBLE_POINT` and `termination_status` is
   neither `INFEASIBLE` nor `INFEASIBLE_OR_UNBOUNDED`. `NEARLY_FEASIBLE_POINT` is accepted
   because stock MES accepts it (`has_values`) and solvers return it for points feasible to
   slightly relaxed tolerances (e.g. barrier without crossover); likewise a feasible point
   with a limit-type termination status is used as before. A warning is logged (at most five
   times) when the accepted status is not `OPTIMAL`/`FEASIBLE_POINT`. Everything else goes to
   phase I, or is an error when `ExpectFeasibleSubproblems = true`.
2. *`ExpectFeasibleSubproblems = true` and no feasible point:* an immediate `error`. MES runs
   `compute_conflict!`, displays the conflict, logs a non-throwing `@error` and then fails with
   an `UndefVarError`; `compute_conflict!` on these models does not return with HiGHS
   (MathOptIIS).
3. *Phase I:* an immediate `error` if the phase-I problem has no (nearly) feasible point, where
   MES tests `has_values`, runs `compute_conflict!` and logs a non-throwing `@error`.

Remove this file once MES is fixed (test "T7 stock MES" in `test/test_nse_budget.jl` then
reports an unexpected pass).
"""
function mes_require_feasible_point!()
    MES_FEASIBLE_POINT_SHIM_APPLIED[] && return true
    return apply_mes_patch!(; force = true)
end

# Called from GenX.__init__ (force = false: honours the escape hatch) and from the public guard.
function apply_mes_patch!(; force::Bool = false)
    MES_FEASIBLE_POINT_SHIM_APPLIED[] && return true
    if !force && get(ENV, MES_PATCH_DISABLE_ENV, "") == "1"
        @info "GenX: $MES_PATCH_DISABLE_ENV=1, MacroEnergySolvers is left unpatched (stock solve_subproblem)."
        return false
    end
    # Evaluating into another module is not allowed while a package that depends on GenX is
    # being precompiled; GenX.__init__ runs again, and patches, when that package is loaded.
    ccall(:jl_generating_output, Cint, ()) == 1 && return false
    v = pkgversion(MacroEnergySolvers)
    if v != MES_PATCHED_VERSION
        if mes_source_has_defect()
            @warn "GenX: MacroEnergySolvers $v is not the version GenX knows how to patch " *
                  "($MES_PATCHED_VERSION), but its solve_subproblem still branches on has_values(m). " *
                  "NOT patched: Benders with infeasible subproblems and HiGHS may return wrong bounds."
        else
            @debug "GenX: MacroEnergySolvers $v: no patch needed or known."
        end
        return false
    end
    @eval MacroEnergySolvers function solve_subproblem(m::Model, planning_sol::NamedTuple,
            linking_variables_sub::Vector{String}, expect_feasible_subproblems::Bool)
        fix_linking_variables!(m, planning_sol, linking_variables_sub)
        optimize!(m)
        if $(mes_has_feasible_point)(m)
            op_cost = objective_value(m)
            lambda = [dual(FixRef(variable_by_name(m, y))) for y in linking_variables_sub]
            theta_coeff = 1
        elseif expect_feasible_subproblems == true
            error("The subproblem did not return a feasible point (termination status " *
                  "$(termination_status(m)), primal status $(primal_status(m))), but " *
                  "ExpectFeasibleSubproblems = true. Set it to false to generate feasibility cuts.")
        else
            @info "Subproblem is infeasible, generating feasibility cut..."
            # Phase-I problem of MES: min slack_max, every affine row relaxed by slack_max
            unfix.(m[:slack_max])
            objfun = objective_function(m)
            @objective(m, Min, m[:slack_max])
            optimize!(m)
            if !$(mes_has_feasible_point)(m)
                error("Feasibility subproblem has no feasible point (termination status " *
                      "$(termination_status(m)), primal status $(primal_status(m))); " *
                      "this should not happen. Check the model.")
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
    @info "GenX: patched MacroEnergySolvers $v solve_subproblem (infeasible-point defect: an " *
          "infeasible subproblem that returns a point was taken for solved). " *
          "Set $MES_PATCH_DISABLE_ENV=1 to keep stock behaviour."
    return true
end

# The status test of the patched method; see `mes_require_feasible_point!`, difference 1.
function mes_has_feasible_point(m::Model)
    ts = termination_status(m)
    ps = primal_status(m)
    (ts == MOI.INFEASIBLE || ts == MOI.INFEASIBLE_OR_UNBOUNDED) && return false
    (ps == MOI.FEASIBLE_POINT || ps == MOI.NEARLY_FEASIBLE_POINT) || return false
    if ts != MOI.OPTIMAL || ps != MOI.FEASIBLE_POINT
        @warn "Benders subproblem accepted with termination status $ts, primal status $ps." maxlog=5
    end
    return true
end

# Does the installed MES still branch on `has_values(m)` right after solving the subproblem?
function mes_source_has_defect()
    f = joinpath(dirname(pathof(MacroEnergySolvers)), "benders", "subproblems.jl")
    isfile(f) || return false
    src = read(f, String)
    i = findfirst("function solve_subproblem(", src)
    isnothing(i) && return false
    return occursin(r"optimize!\(m\)\s*if has_values\(m\)", src[first(i):end])
end
