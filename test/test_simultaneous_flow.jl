module TestSimultaneousFlow

# Tests of the optional directional line ratings (`Profile_Forward`, `Profile_Reverse` columns of
# Network.csv) and the optional simultaneous-flow groups (system/Simultaneous_Flow_Constraints.csv),
# ported from the chair's GenX fork (/home/gm1710/GenX @ 49557d5c, transmission.jl). Case:
# test/benders/1_three_zones (3 zones, lines 1 = MA->CT and 2 = MA->ME, both loss lines,
# Trans_Loss_Segments 1, ParameterScale 1, 4 representative periods x 6 h). Every variant is a
# COPY of the case in a temporary folder with the new file/columns written there; the repo case
# is never written to.
#
# Model (lines l, hours t, loss-line flow split vFLOW = vTAUX_POS - vTAUX_NEG):
#   vTAUX_POS[l,t] <= eAvail_Trans_Cap[l] * Profile_Forward[l]      (only if the columns exist)
#   vTAUX_NEG[l,t] <= eAvail_Trans_Cap[l] * Profile_Reverse[l]      (only if the columns exist)
#   sum_{(l,d) in group g} (d == forward ? vTAUX_POS[l,t] : vTAUX_NEG[l,t]) <= limit_g  (per g, t)

include(joinpath(@__DIR__, "nse_budget_helpers.jl"))

const MW = 1.0  # limits in the file are MW; the model is in GW when ParameterScale = 1

"Copy of the case with optional profile columns (Vector of (fwd, rev) per line) and simflow file text."
function sf_case(; profiles = nothing, simflow = nothing, loss_override = nothing,
        overrides = Dict())
    dir = joinpath(mktempdir(), "case")
    cp(CASE, dir)
    for d in filter(startswith("results"), readdir(dir))
        rm(joinpath(dir, d); recursive = true)
    end
    netf = joinpath(dir, "system", "Network.csv")
    lines = readlines(netf)
    if loss_override !== nothing   # (line number, new Line_Loss_Percentage string)
        l, val = loss_override
        cells = split(lines[1 + l], ",")
        hdr = split(lines[1], ",")
        cells[findfirst(==("Line_Loss_Percentage"), hdr)] = val
        lines[1 + l] = join(cells, ",")
    end
    if profiles !== nothing
        lines[1] *= ",Profile_Forward,Profile_Reverse"
        for i in 2:length(lines)
            l = i - 1
            lines[i] *= l <= length(profiles) ? ",$(profiles[l][1]),$(profiles[l][2])" : ",,"
        end
    end
    write(netf, join(lines, "\n") * "\n")
    if simflow !== nothing
        write(joinpath(dir, "system", "Simultaneous_Flow_Constraints.csv"),
            "Simultaneous Flow Group,Transmission Path,Line_Number,Direction,limit_MW\n" * simflow)
    end
    if !isempty(overrides)
        f = joinpath(dir, "settings", "genx_settings.yml")
        y = YAML.load_file(f)
        merge!(y, overrides)
        YAML.write_file(f, y)
    end
    return dir
end

load_case(dir, setup) = quiet(() -> GenX.load_inputs(setup, dir))

function directional(EP)
    pos = value.(EP[:vTAUX_POS])
    neg = value.(EP[:vTAUX_NEG])
    T = size(value.(EP[:vFLOW]), 2)
    P = [pos[l, t] for l in 1:2, t in 1:T]
    N = [neg[l, t] for l in 1:2, t in 1:T]
    return P, N
end

cname_count(EP, prefix) = count(c -> startswith(name(c), prefix),
    all_constraints(EP; include_variable_in_set_constraints = false))

"The row named `nm` as built: (number of terms, [coefficient of each variable in `vars`], rhs); printed."
function row_terms(EP, nm, vars)
    c = constraint_by_name(EP, nm)
    c === nothing && return (-1, fill(NaN, length(vars)), NaN)
    f = constraint_object(c).func
    println("row $nm : $f <= $(normalized_rhs(c))")
    return (length(f.terms), [normalized_coefficient(c, v) for v in vars], normalized_rhs(c))
end

"The message of the error thrown by `f()`, or \"\" if it throws none."
errmsg(f) =
    try
        f()
        ""
    catch e
        sprint(showerror, e)
    end

const SETUP = case_setup()
const SC = scale(SETUP)
const BASE_DIR = sf_case()
const INPUTS_BASE = load_case(BASE_DIR, SETUP)
const EP_BASE = solve_monolithic(SETUP, INPUTS_BASE)
const Z_BASE = objective_value(EP_BASE)
const P_BASE, N_BASE = directional(EP_BASE)
println("SimFlow tests: base objective $Z_BASE; max fwd flow per line (MW) ",
    [maximum(P_BASE[l, :]) * SC for l in 1:2], "; max rev flow per line (MW) ",
    [maximum(N_BASE[l, :]) * SC for l in 1:2])

@testset "Simultaneous flow and directional ratings" begin
    @testset "S0 absent file and columns: today's model" begin
        @test Z_BASE≈KNOWN_OPTIMUM_CASE1 rtol=OBJECTIVE_RTOL
        @test !haskey(object_dictionary(EP_BASE), :cTAuxPosProfile)
        @test !haskey(object_dictionary(EP_BASE), :cTAuxNegProfile)
        @test cname_count(EP_BASE, "cSimFlow") == 0
        # the unmodified repo case gives the identical model and objective
        inputs_repo = case_inputs(SETUP; single_segment = false)
        EP_repo = solve_monolithic(SETUP, inputs_repo)
        @test num_constraints(EP_repo; count_variable_in_set_constraints = true) ==
              num_constraints(EP_BASE; count_variable_in_set_constraints = true)
        @test num_variables(EP_repo) == num_variables(EP_BASE)
        @test objective_value(EP_repo) == Z_BASE
        # all-ones profiles and a slack group: redundant rows, same optimum
        dir = sf_case(profiles = [(1.0, 1.0), (1.0, 1.0)],
            simflow = "Slack,MA_to_CT,1,forward,1e6\nSlack,MA_to_ME,2,forward,1e6\n")
        EP1 = solve_monolithic(SETUP, load_case(dir, SETUP))
        @test haskey(object_dictionary(EP1), :cTAuxPosProfile)
        @test cname_count(EP1, "cSimFlowSlack") == INPUTS_BASE["T"]
        @test objective_value(EP1)≈Z_BASE rtol=LP_RTOL
    end

    @testset "S1 a group limit below the unconstrained sum binds" begin
        # group = both lines of zone MA in the direction that carries the larger peak sum
        sf = vec(sum(P_BASE; dims = 1))
        sr = vec(sum(N_BASE; dims = 1))
        dir_name, s = maximum(sf) >= maximum(sr) ? ("forward", sf) : ("reverse", sr)
        @test maximum(s) > 0
        limit_mw = 0.5 * maximum(s) * SC
        println("S1: group = lines 1,2 $dir_name; unconstrained peak sum $(maximum(s) * SC) MW; limit $limit_mw MW")
        dir = sf_case(simflow = "G,MA_to_CT,1,$dir_name,$limit_mw\nG,MA_to_ME,2,$dir_name,$limit_mw\n")
        inputs = load_case(dir, SETUP)
        EP = solve_monolithic(SETUP, inputs)
        P, N = directional(EP)
        g = vec(sum(dir_name == "forward" ? P : N; dims = 1)) .* SC
        tol = 1e-6 * limit_mw
        @test all(g .<= limit_mw + tol)
        nbind = count(x -> x >= limit_mw - tol, g)
        println("S1: hours at the limit $nbind of $(length(g)); objective $(objective_value(EP)) vs base $Z_BASE")
        @test nbind >= 1
        @test objective_value(EP) > Z_BASE * (1 + 1e-7)
        @test cname_count(EP, "cSimFlowG") == inputs["T"]
    end

    @testset "S2 a one-way profile forbids the other direction" begin
        # the (line, direction) with the largest unconstrained flow is closed
        cand = [(maximum(P_BASE[l, :]), l, :fwd) for l in 1:2]
        append!(cand, [(maximum(N_BASE[l, :]), l, :rev) for l in 1:2])
        fmax, l0, d0 = maximum(cand)
        @test fmax > 0
        prof = [(1.0, 1.0), (1.0, 1.0)]
        prof[l0] = d0 == :fwd ? (0.0, 1.0) : (1.0, 0.0)
        println("S2: closing line $l0 direction $d0 (unconstrained peak $(fmax * SC) MW)")
        EP = solve_monolithic(SETUP, load_case(sf_case(profiles = prof), SETUP))
        P, N = directional(EP)
        closed = d0 == :fwd ? P[l0, :] : N[l0, :]
        flow = value.(EP[:vFLOW])[l0, :]
        @test maximum(closed) * SC <= 1e-6
        @test (d0 == :fwd ? maximum(flow) : -minimum(flow)) * SC <= 1e-6
        @test objective_value(EP) >= Z_BASE * (1 - LP_RTOL)
        # a partial rating 0.5 caps the open direction at half the line
        prof2 = [(1.0, 1.0), (1.0, 1.0)]
        prof2[l0] = d0 == :fwd ? (0.5, 1.0) : (1.0, 0.5)
        inputs2 = load_case(sf_case(profiles = prof2), SETUP)
        EP2 = solve_monolithic(SETUP, inputs2)
        P2, N2 = directional(EP2)
        cap = value(EP2[:eAvail_Trans_Cap][l0])
        @test maximum(d0 == :fwd ? P2[l0, :] : N2[l0, :]) <= 0.5 * cap * (1 + 1e-6) + 1e-9
    end

    @testset "S3 refusals" begin
        function load_and_build(dir)
            inputs = load_case(dir, SETUP)
            s = copy(SETUP)
            s["Benders"] = 0
            quiet() do
                GenX.generate_model(s, inputs, GenX.configure_solver(SETTINGS, HiGHS.Optimizer))
            end
        end
        bad(simflow; kw...) = @test_throws Exception load_and_build(sf_case(; simflow, kw...))
        bad("G,MA_to_CT,1,both,100\n")                                   # direction
        bad("G,MA_to_CT,3,forward,100\n")                                # line > L
        bad("G,MA_to_CT,0,forward,100\n")                                # line < 1
        bad("G,MA_to_CT,1,forward,100\nG,MA_to_ME,2,forward,200\n")      # two limits in one group
        bad("G,MA_to_ME,2,forward,100\n"; loss_override = (2, "0"))      # not a loss line
        # the same (line, direction) twice in one group: the loader's message
        msg = errmsg(() -> load_case(
            sf_case(simflow = "G,MA_to_CT,1,forward,100\nG,MA_to_CT,1,forward,100\n"), SETUP))
        @test occursin("Simultaneous_Flow_Constraints.csv: group G lists line 1 forward twice", msg)
        # Trans_Loss_Segments other than 1 (0 here: no Ohms/kV columns needed) with a profile
        # column or with a group file
        s0 = copy(SETUP)
        s0["Trans_Loss_Segments"] = 0
        msg = errmsg(() -> load_case(sf_case(profiles = [(1.0, 1.0), (1.0, 1.0)]), s0))
        @test occursin(
            "Network.csv: Profile_Forward/Profile_Reverse are implemented for Trans_Loss_Segments = 1 only",
            msg)
        msg = errmsg(() -> load_case(sf_case(simflow = "G,MA_to_CT,1,forward,100\n"), s0))
        @test occursin(
            "Simultaneous_Flow_Constraints.csv: simultaneous-flow groups are implemented for Trans_Loss_Segments = 1 only",
            msg)
    end

    @testset "S4 Benders: rows in the subproblems only, bounds bracket the monolith" begin
        sf = vec(sum(P_BASE; dims = 1))
        sr = vec(sum(N_BASE; dims = 1))
        dir_name, s = maximum(sf) >= maximum(sr) ? ("forward", sf) : ("reverse", sr)
        limit_mw = 0.5 * maximum(s) * SC
        txt = "G,MA_to_CT,1,$dir_name,$limit_mw\nG,MA_to_ME,2,$dir_name,$limit_mw\n"
        dir = sf_case(simflow = txt, profiles = [(1.0, 0.75), (1.0, 1.0)])
        inputs = load_case(dir, SETUP)
        z_mono = objective_value(solve_monolithic(SETUP, inputs))
        bi = build_benders(benders_setup(SETUP), inputs)
        @test cname_count(bi["planning_problem"], "cSimFlow") == 0
        @test !haskey(object_dictionary(bi["planning_problem"]), :cTAuxNegProfile)
        nsub = [cname_count(sp[:model], "cSimFlowG") for sp in bi["subproblems"]]
        @test sum(nsub) == inputs["T"]
        @test all(n -> n == inputs["hours_per_subperiod"], nsub)
        @test all(sp -> haskey(object_dictionary(sp[:model]), :cTAuxNegProfile), bi["subproblems"])
        # end-to-end Benders on the same files
        bdir = sf_case(simflow = txt, profiles = [(1.0, 0.75), (1.0, 1.0)],
            overrides = Dict("Benders" => 1, "OverwriteResults" => 1, "PrintModel" => 0))
        redirect_stdout(devnull) do
            with_logger(ConsoleLogger(stderr, Logging.Warn)) do
                run_genx_case!(bdir, HiGHS.Optimizer)
            end
        end
        conv = CSV.read(joinpath(bdir, "results_benders", "benders_convergence.csv"), DataFrame)
        ub, lb = conv[end, :UB], conv[end, :LB]
        println("S4: monolith $z_mono; Benders LB $lb UB $ub; base $Z_BASE")
        @test z_mono > Z_BASE * (1 + 1e-7)
        @test lb <= z_mono * (1 + 1e-6)
        @test ub >= z_mono * (1 - 1e-6)
        @test abs(ub - z_mono) / abs(z_mono) <= OBJECTIVE_RTOL
    end

    # S5 and S6 use FIXED (hard-coded) directions. Line 2 = MA_to_ME, Start_Zone 1, End_Zone 3:
    # its forward component vTAUX_POS[2,t] is flow from zone 1 to zone 3. Line 1 = MA_to_CT,
    # zones 1 -> 2: its reverse component vTAUX_NEG[1,t] is flow from zone 2 to zone 1.
    # Unconstrained solution of this case (printed at the top; della jobs 18302937, 18273509):
    # peak forward flow on line 2 = 1035.7 MW, peak reverse flow on line 1 = 5900 MW
    # (= 2950 existing + 2950 maximum reinforcement), peak forward flow on line 1 = 0.
    # The tests assert the parts of this that they rely on.

    @testset "S5 forward directional rating (line 2): row as built, binds, 0 closes" begin
        @test maximum(P_BASE[2, :]) * SC > 1000.0
        prof = [(0.5, 1.0), (0.25, 1.0)]
        inputs = load_case(sf_case(profiles = prof), SETUP)
        EP = solve_monolithic(SETUP, inputs)
        T = inputs["T"]
        # (i) rows by name. ParameterScale 1: the model is in GW. eAvail_Trans_Cap[l] =
        # Line_Max_Flow_MW/1000 + vNEW_TRANS_CAP[l]; the row vTAUX_POS <= Profile_Forward * eAvail
        # is stored as vTAUX_POS - pf * vNEW_TRANS_CAP <= pf * existing:
        #   line 1: 2950 MW = 2.95 GW, pf 0.5:  vTAUX_POS[1,1] - 0.5 vNEW_TRANS_CAP[1]  <= 0.5 * 2.95 = 1.475
        #   line 2: 2000 MW = 2.0 GW,  pf 0.25: vTAUX_POS[2,1] - 0.25 vNEW_TRANS_CAP[2] <= 0.25 * 2.0 = 0.5
        # and the reverse rows with pr 1.0: vTAUX_NEG[2,1] - vNEW_TRANS_CAP[2] <= 2.0
        n, a, rhs = row_terms(EP, "cTAuxPosProfile[1,1]",
            [EP[:vTAUX_POS][1, 1], EP[:vNEW_TRANS_CAP][1]])
        @test n == 2
        @test a == [1.0, -0.5]
        @test rhs ≈ 1.475 atol = 1e-12
        n, a, rhs = row_terms(EP, "cTAuxPosProfile[2,1]",
            [EP[:vTAUX_POS][2, 1], EP[:vNEW_TRANS_CAP][2]])
        @test n == 2
        @test a == [1.0, -0.25]
        @test rhs ≈ 0.5 atol = 1e-12
        n, a, rhs = row_terms(EP, "cTAuxNegProfile[2,1]",
            [EP[:vTAUX_NEG][2, 1], EP[:vNEW_TRANS_CAP][2]])
        @test n == 2
        @test a == [1.0, -1.0]
        @test rhs ≈ 2.0 atol = 1e-12
        @test cname_count(EP, "cTAuxPosProfile") == 2 * T
        # (ii) the rating binds. Line 2's available capacity is at most 2000 + 2000 (maximum
        # reinforcement) = 4000 MW, so its forward cap is at most 0.25 * 4000 = 1000 MW, below the
        # unconstrained peak of 1035.7 MW: the cap cuts the unconstrained flow, and the forward
        # flow sits at 0.25 * capacity in at least one hour.
        P, N = directional(EP)
        cap2 = value(EP[:eAvail_Trans_Cap][2])
        println("S5: line 2 capacity $(cap2 * SC) MW; forward cap $(0.25 * cap2 * SC) MW; peak forward ",
            "$(maximum(P[2, :]) * SC) MW (unconstrained $(maximum(P_BASE[2, :]) * SC) MW); objective ",
            "$(objective_value(EP)) vs base $Z_BASE")
        @test cap2 <= 4.0 + 1e-9
        @test maximum(P_BASE[2, :]) > 0.25 * cap2 * (1 + 1e-6)
        @test all(P[2, :] .<= 0.25 * cap2 * (1 + 1e-7) + 1e-9)
        @test maximum(P[2, :])≈0.25 * cap2 rtol=1e-6
        @test objective_value(EP) > Z_BASE * (1 + 1e-7)
        # (iii) a forward rating of 0 on line 2: no forward flow on line 2 in any hour
        EP0 = solve_monolithic(SETUP, load_case(sf_case(profiles = [(1.0, 1.0), (0.0, 1.0)]), SETUP))
        P0, N0 = directional(EP0)
        @test maximum(P0[2, :]) * SC <= 1e-6
        @test maximum(value.(EP0[:vFLOW])[2, :]) * SC <= 1e-6
        @test objective_value(EP0) > Z_BASE * (1 + 1e-7)
    end

    @testset "S6 hard-coded forward group and mixed group: rows, binding, terms add" begin
        T = INPUTS_BASE["T"]
        # Group Fwd: one forward term, line 2, 600 MW. Group Mixed: line 2 forward + line 1
        # reverse, 3000 MW. ParameterScale 1: limits 600/1000 = 0.6 GW and 3000/1000 = 3.0 GW.
        #   cSimFlowFwd[t]:   vTAUX_POS[2,t]                  <= 0.6
        #   cSimFlowMixed[t]: vTAUX_POS[2,t] + vTAUX_NEG[1,t] <= 3.0
        # (i) both rows as built, from one case holding both groups
        both = "Fwd,MA_to_ME,2,forward,600\nMixed,MA_to_ME,2,forward,3000\nMixed,MA_to_CT,1,reverse,3000\n"
        inputs = load_case(sf_case(simflow = both), SETUP)
        EPb = solve_monolithic(SETUP, inputs)
        @test cname_count(EPb, "cSimFlowFwd") == T
        @test cname_count(EPb, "cSimFlowMixed") == T
        for t in (1, T)
            n, a, rhs = row_terms(EPb, "cSimFlowFwd[$t]", [EPb[:vTAUX_POS][2, t]])
            @test n == 1
            @test a == [1.0]
            @test rhs ≈ 0.6 atol = 1e-12
            n, a, rhs = row_terms(EPb, "cSimFlowMixed[$t]",
                [EPb[:vTAUX_POS][2, t], EPb[:vTAUX_NEG][1, t]])
            @test n == 2
            @test a == [1.0, 1.0]
            @test rhs ≈ 3.0 atol = 1e-12
        end

        # (ii) the forward limit binds (group Fwd alone): unconstrained peak 1035.7 MW > 600 MW
        @test maximum(P_BASE[2, :]) * SC > 600
        EPf = solve_monolithic(SETUP, load_case(sf_case(simflow = "Fwd,MA_to_ME,2,forward,600\n"), SETUP))
        Pf, Nf = directional(EPf)
        println("S6 Fwd: peak forward flow on line 2 $(maximum(Pf[2, :]) * SC) MW (limit 600; unconstrained ",
            "$(maximum(P_BASE[2, :]) * SC)); hours at the limit $(count(x -> x * SC >= 600 - 1e-4, Pf[2, :])) of $T")
        @test all(Pf[2, :] .* SC .<= 600 + 1e-4)
        @test maximum(Pf[2, :]) * SC≈600 rtol=1e-6
        @test objective_value(EPf) > Z_BASE * (1 + 1e-7)

        # (iii) group Mixed alone: the two gross components ADD.
        # Unconstrained: peak line-1 reverse 5900 MW, peak line-2 forward 1035.7 MW.
        # Gross form (as built): P2 + N1 <= 3000 in every hour, so N1 <= 3000 - P2 <= 3000 MW:
        #   the line-1 reverse peak must fall from 5900 MW to at most 3000 MW.
        # A net form (P2 - N1 <= 3000, a reverse flow offsetting the forward flow) would be slack
        #   at the unconstrained solution in every hour (P2 - N1 <= P2 <= 1035.7 < 3000), leaving
        #   the line-1 reverse peak at 5900 MW and the objective at the base value.
        @test maximum(P_BASE[2, :]) * SC < 3000
        @test maximum(N_BASE[1, :]) * SC > 3000
        EPm = solve_monolithic(SETUP,
            load_case(sf_case(simflow = "Mixed,MA_to_ME,2,forward,3000\nMixed,MA_to_CT,1,reverse,3000\n"), SETUP))
        Pm, Nm = directional(EPm)
        gsum = (Pm[2, :] .+ Nm[1, :]) .* SC
        tb = argmax(gsum)
        println("S6 Mixed: peak gross sum $(gsum[tb]) MW at hour $tb (line-2 forward $(Pm[2, tb] * SC) + line-1 ",
            "reverse $(Nm[1, tb] * SC)); peak line-1 reverse $(maximum(Nm[1, :]) * SC) MW; hours with both ",
            "components > 1 MW: $(count(t -> Pm[2, t] * SC > 1 && Nm[1, t] * SC > 1, 1:T)); objective ",
            "$(objective_value(EPm)) vs base $Z_BASE")
        @test all(gsum .<= 3000 + 1e-4)
        @test gsum[tb]≈3000 rtol=1e-6
        @test maximum(Nm[1, :]) * SC <= 3000 + 1e-4
        @test objective_value(EPm) > Z_BASE * (1 + 1e-7)
    end
end

end # module TestSimultaneousFlow
