using Base.Test

println("test MUMPS")
@testset "MUMPS" begin
    @testset "DivGrad" begin include("testDivGrad.jl") end
    @testset "Two systems" begin include("testTwoSystem.jl") end
    @testset "Two systems parallel" begin include("testTwoSystemParallel.jl") end
    @testset "Sparse RHS" begin include("testSparseRHS.jl") end
    println("Done!")
end
