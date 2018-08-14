
using Test
using Distributed
addprocs(1)

# Factorize two different matrices (one real, one complex) and solve for multiple rhs and free memory
@everywhere using MUMPS

include("getDivGrad.jl")

Ar = getDivGrad(30, 31, 33)
n  = size(Ar, 1)
Ac  = Ar + im * spdiagm(0 => rand(n))

nrhs = 10
rhsr = randn(n, nrhs) + im * randn(n, nrhs)
rhsc = randn(n, nrhs)

println("Factorize real matrix on main instance and destroy")
fact = factorMUMPS(Ar, 1)
destroyMUMPS(fact)

println("Factorize real matrix on remote worker and destroy")
fact = remotecall_fetch(factorMUMPS, workers()[1], Ar, 1)
remotecall_fetch(destroyMUMPS, workers()[1], fact)
