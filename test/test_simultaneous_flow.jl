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
end

end # module TestSimultaneousFlow
