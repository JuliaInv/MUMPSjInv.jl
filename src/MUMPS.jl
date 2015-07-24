module MUMPS

type MUMPSfactorization{T}
	ptr::Int64     # pointer to factorization
	worker::Int64  # id of worker that holds factorization
	n::Int64       # matrix size
	a11::T         # first element (HACK for generating the right parametric type)
	time::Float64  # factorization time
end
	
	include("MUMPSfuncs.jl")
	
	export solveMUMPS,solveMUMPS!, factorMUMPS, applyMUMPS,applyMUMPS!,destroyMUMPS, MUMPSfactorization
	
end
