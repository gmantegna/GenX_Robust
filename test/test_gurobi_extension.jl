module TestGurobiExtension

using Test
using GenX

# Regression test: GenX itself must NOT own a method of `benders_gurobi_optimizer`.
# A stub method in GenX with the same signature as the one in `GenXGurobiExt` makes the
# extension overwrite it, and Julia (>= 1.10) then refuses to precompile the extension
# ("Method overwriting is not permitted during Module precompilation").
@testset "benders_gurobi_optimizer is only implemented by GenXGurobiExt" begin
    ms = methods(GenX.benders_gurobi_optimizer)
    @test all(m -> m.module !== GenX, ms)
    @test length(ms) <= 1

    if isempty(ms)
        # Gurobi not loaded: no method, and the MethodError carries an informative hint.
        @test isnothing(Base.get_extension(GenX, :GenXGurobiExt))
        @test_throws MethodError GenX.benders_gurobi_optimizer(Dict())
        msg = try
            GenX.benders_gurobi_optimizer(Dict())
        catch e
            sprint(showerror, e)
        end
        @test occursin("using Gurobi", msg)
    else
        # Compare by name: `Base.get_extension` returns `nothing` if the extension's
        # `__init__` failed (e.g. no Gurobi license) even though its methods are defined.
        @test nameof(only(ms).module) === :GenXGurobiExt
    end
end

end # module TestGurobiExtension
