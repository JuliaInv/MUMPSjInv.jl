using Test

@testset "MUMPSjInv" begin
    @testset "DivGrad" begin include("testDivGrad.jl") end
    @testset "Two systems" begin include("testTwoSystem.jl") end
    @testset "Two systems parallel" begin include("testTwoSystemParallel.jl") end
    @testset "Sparse RHS" begin include("testSparseRHS.jl") end
    @testset "testDestroyMUMPS" begin include("testDestroyMUMPS.jl") end
end
