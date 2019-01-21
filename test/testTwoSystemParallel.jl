# Factorize two different matrices (one real, one complex) and solve for multiple rhs and free memory
using Distributed
if nworkers()<3
    addprocs(3)
end
@everywhere using MUMPS
using Test

include("getDivGrad.jl");
A = getDivGrad(24,20,22);
A2 = getDivGrad(30,31,33);
n = size(A,1);
n2 = size(A2,1);
A  = A + im*spdiagm(0 => rand(n));
nrhs = 10;
rhs = randn(n,nrhs) + im*randn(n,nrhs);
rhs2 = randn(n2,nrhs);

@test nworkers() >= 3
println("Factorize complex matrix on worker 1")
@time F1 = remotecall_fetch(factorMUMPS,workers()[1], A, 1)
dump(F1)

println("Factorize real matrix on worker 2")
@time F2 = remotecall_fetch(factorMUMPS, workers()[2], A2, 1)
dump(F2)

println("Solve complex system on worker 1")
@time x = remotecall_fetch(applyMUMPS, workers()[3], F1, rhs)
err = zeros(nrhs)
for i=1:nrhs
    err[i] =  norm(A * x[:, i]-rhs[:, i]) / norm(rhs[:, i])
end
@test maximum(err) < 1e-14

println("Solve real system on worker 2")
@time x2 = remotecall_fetch(applyMUMPS, workers()[2], F2, rhs2)
err = zeros(nrhs)
for i = 1:nrhs
    err[i] = norm(A2 * x2[:, i] - rhs2[:, i]) / norm(rhs2[:, i])
end
@test maximum(err) < 1e-14

println("Free memory on worker 1")
destroyMUMPS(F1)

println("Free memory on worker 2")
destroyMUMPS(F2)
