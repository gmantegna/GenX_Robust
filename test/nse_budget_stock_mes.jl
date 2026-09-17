# Run by test_nse_budget.jl in a SUBPROCESS with GENX_DISABLE_MES_PATCH=1, i.e. with stock
# MacroEnergySolvers. Not a test file by itself: it prints KEY=value lines that the parent
# parses. Two questions:
#   (a) does stock MES still take an infeasible subproblem for a solved one (HiGHS)? If the
#       answer becomes "no", MES has been fixed upstream and GenX's patch can be deleted;
#   (b) does run_genx_case! refuse to write outputs when Benders ends with a negative gap?
include(joinpath(@__DIR__, "nse_budget_helpers.jl"))

println("PATCH_APPLIED=", isdefined(GenX, :MES_FEASIBLE_POINT_SHIM_APPLIED) ?
                          GenX.MES_FEASIBLE_POINT_SHIM_APPLIED[] : "undefined")

function stock_benders()
    s0 = case_setup()
    inputs = case_inputs(s0)
    estar = expected_shed_mwh(solve_monolithic(s0, inputs), inputs, s0)
    s = budget_setup(0.5 * estar)
    z_mono = objective_value(solve_monolithic(s, inputs))
    results, nfeas, _, _, _ = solve_benders(s, inputs; resolve_incumbent = false)
    println("STOCK_ZMONO=", z_mono)
    println("STOCK_STATUS=", results.termination_status)
    println("STOCK_UB=", results.UB_hist[end])
    println("STOCK_LB=", results.LB_hist[end])
    println("STOCK_ITERS=", length(results.UB_hist))
    println("STOCK_FEASCUTS=", nfeas)
end
stock_benders()

function stock_runner()
    case = e2e_benders_case()
    msg = "none"
    try
        redirect_stdout(devnull) do
            with_logger(ConsoleLogger(stderr, Logging.Warn)) do
                GenX.run_genx_case!(case, HiGHS.Optimizer)
            end
        end
    catch e
        msg = replace(sprint(showerror, e), '\n' => ' ')
        msg = string(typeof(e), ": ", first(msg, 300))
    end
    println("RUNNER_ERROR=", msg)
    println("RUNNER_RESULTS_WRITTEN=", any(startswith("results"), readdir(case)))
end
stock_runner()
println("STOCK DONE")
