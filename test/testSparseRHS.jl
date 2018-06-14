using MUMPS
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
include("getDivGrad.jl");

A = getDivGrad(32,32,16);
n = size(A,1);

# REAL: test for single rhs
println("Test for real SPD matrix: one sparse random rhs");
rhs = sprandn(n,1,0.05);

x = solveMUMPS(A,rhs,1);
err  =  norm(A*x-full(rhs)) / norm(full(rhs));
@test err < 1e-14

# rhs as sparse vector instead of n X 1 sparse matrix
rhs = vec(rhs)

x = solveMUMPS(A,rhs,1);
err  =  norm(A*x-full(rhs)) / norm(full(rhs));
@test err < 1e-14

# REAL: test for multiple rhs
println("Test for real SPD matrix: multiple sparse random rhs");
nrhs = 10;
rhs = sprandn(n,nrhs,0.03);

x = solveMUMPS(A,rhs,1);

err = zeros(nrhs)
for i=1:nrhs
        err[i] =  norm(A*x[:,i]-full(rhs[:,i])) / norm(full(rhs[:,i]));
end
@test eltype(x) == Float64
@test maximum(err) < 1e-14

println("Test for real SPD matrix: multiple sparse delta function rhs");
rhs = speye(n)
rhs = rhs[:,1:nrhs]

x = solveMUMPS(A,rhs,1);

err = zeros(nrhs)
for i=1:nrhs
        err[i] =  norm(A*x[:,i]-full(rhs[:,i])) / norm(full(rhs[:,i]));
end
@test eltype(x) == Float64
@test maximum(err) < 1e-14

# COMPLEX: test for one rhs
println("Test for complex SPD matrix: one rhs");
r = rand(n);
A = A + im*spdiagm(r,0);

rhs = sprandn(n,1,0.03) + im*sprandn(n,1,0.03)

x = solveMUMPS(A,rhs,2);

err=  norm(A*x-full(rhs)) / norm(full(rhs))
@test eltype(x) == Complex128
@test err < 1e-14

# sparse vector rhs
rhs = vec(sprandn(n,1,0.03) + im*sprandn(n,1,0.03))

x = solveMUMPS(A,rhs,2);

err=  norm(A*x-full(rhs)) / norm(full(rhs))
@test eltype(x) == Complex128
@test err < 1e-14

# COMPLEX : test for multiple rhs
println("Test for complex SPD matrix: multiple sparse random rhs");
nrhs = 10;
rhs = sprandn(n,nrhs,0.03) + im*sprandn(n,nrhs,0.03);

x = solveMUMPS(A,rhs,2);
err = zeros(nrhs)
for i=1:nrhs
        err[i] =  norm(A*x[:,i]-full(rhs[:,i])) / norm(full(rhs[:,i]));
end
@test eltype(x) == Complex128
@test maximum(err) < 1e-14

println("DONE!")
