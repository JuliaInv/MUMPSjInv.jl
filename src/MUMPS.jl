module MUMPS
	abstract MUMPSfactorization
	type MUMPSfactorizationReal <: MUMPSfactorization
		ptr::Int64     # pointer to factorization
		worker::Int64  # id of worker that holds factorization
		n::Int64       # matrix size
		time::Float64  # factorization time
	end

	type MUMPSfactorizationComplex <: MUMPSfactorization
		ptr::Int64     # pointer to factorization
		worker::Int64  # id of worker that holds factorization
		n::Int64       # matrix size
		time::Float64  # factorization time
	end
	
	arrayOrSparseCSC        = Union{Array,SparseMatrixCSC}
	arrayOrSparseCSCReal    = Union{Array{Float64},SparseMatrixCSC{Float64,Int64}}
	arrayOrSparseCSCComplex = Union{Array{Complex{Float64}},SparseMatrixCSC{Complex{Float64},Int64}}
	
	include("MUMPSfuncs.jl")
	
	export solveMUMPS, factorMUMPS, applyMUMPS,destroyMUMPS, MUMPSfactorizationReal,MUMPSfactorizationComplex
	
end
