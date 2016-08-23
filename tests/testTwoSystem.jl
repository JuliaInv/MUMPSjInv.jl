# Factorize two different matrices (one real, one complex) and solve for multiple rhs and free memory

using MUMPS
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
include("getDivGrad.jl");


A  = getDivGrad(24,23,25);
A2 = getDivGrad(34,32,36);

n  = size(A,1);
n2 = size(A2,1);
A  = A + im*spdiagm(rand(n),0);

nrhs = 10;
rhs = randn(n,nrhs) + im*randn(n,nrhs);
rhs2 = randn(n2,nrhs);

println("Factorize complex matrix")
tic();
F1 = factorMUMPS(A,1);
toc();

println("Factorize real matrix")
tic();
F2 = factorMUMPS(A2,1);
toc();

println("Solve complex system")
tic();
x = applyMUMPS(F1,rhs);
toc();
err = zeros(nrhs)
for i=1:nrhs
        err[i] =  norm(A*x[:,i]-rhs[:,i]) / norm(rhs[:,i]);
end
@test maximum(err) < 1e-14

println("Solve real system")
tic();
x2 = applyMUMPS(F2,rhs2);
toc();
err = zeros(nrhs)
for i=1:nrhs
        err[i] =  norm(A2*x2[:,i]-rhs2[:,i]) / norm(rhs2[:,i]);
end
@test maximum(err) < 1e-14

println("Free memory")
destroyMUMPS(F1);
destroyMUMPS(F2);

println("DONE!")


