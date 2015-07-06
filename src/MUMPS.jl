module MUMPS
	
	type MUMPSfactorization
		ptr::Int64     # pointer to factorization
		worker::Int64  # id of worker that holds factorization
		n::Int64       # matrix size
		real::Bool     # whether matrix is real or not
		time::Float64  # factorization time
	end
	
	arrayOrSparseCSC = Union{Array,SparseMatrixCSC}
	
	include("MUMPSfuncs.jl")
	
	export solveMUMPS, factorMUMPS, applyMUMPS,destroyMUMPS, MUMPSfactorization
	
end
