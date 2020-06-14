using DG_Playground, Test, LinearAlgebra


# Check that the derivative of a linear term is constant
@testset "Simple Derivative Test" begin
    for n in 1:7
        α = β = 0.0
        # compute Gauss Lobatto grid
        r = jacobiGL(α, β, n)
        # build differentiation matrix
        D = dmatrix(r, α, β, n)
        bool = norm(D * r .- 1.0) < eps(1.0) * 100
        @test bool
    end
    # Note that n = 8 will fail due to ill-conditioning
end
