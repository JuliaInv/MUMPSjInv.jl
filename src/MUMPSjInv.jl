__precompile__()

module MUMPSjInv

using jInv
using LinearAlgebra
using SparseArrays
using Distributed
using Printf
using Libdl
import jInv.LinearSolvers


mutable struct MUMPSfactorization{T}
	ptr::Int64     # pointer to factorization
	worker::Int64  # id of worker that holds factorization
	n::Int64       # matrix size
	a11::T         # first element (HACK for generating the right parametric type)
	time::Float64  # factorization time
end
	const MUMPSlibPath  = abspath(joinpath(splitdir(Base.source_path())[1],"..","lib","MUMPS"))
	const hasMUMPS = find_library([MUMPSlibPath])!=""
	println("MUMPSlibPath = ",MUMPSlibPath)
	println("hasMUMPS = ",hasMUMPS)

	include("MUMPSfuncs.jl")
	include("mumpsWrapper.jl")

	export solveMUMPS,solveMUMPS!, factorMUMPS, applyMUMPS,applyMUMPS!,destroyMUMPS, MUMPSfactorization

end
