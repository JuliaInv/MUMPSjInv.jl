
function solveMUMPS{T1,T2}(A::SparseMatrixCSC{T1}, rhs::AbstractArray{T2}, sym::Int=0,ooc::Int=0,tr::Int=0)

	x = zeros(promote_type(T1,T2),size(rhs))

	return solveMUMPS!(A,rhs,x,sym,ooc,tr)
end

function solveMUMPS!(A::SparseMatrixCSC, rhs::AbstractArray, x::Array, sym::Int=0,ooc::Int=0,tr::Int=0)

# factorization
factor = factorMUMPS(A,sym,ooc)

# solve system
x = applyMUMPS!(factor,rhs,x,tr)

# free memory
destroyMUMPS(factor)
return x
end

function factorMUMPS(A::SparseMatrixCSC{Complex128},sym::Int=0,ooc::Int=0)
    # Generate LU-factorization of complex matrix A, 'sym' can be: 0=unsymmetric, 1=symm. pos def, 2=general symmetric
    if size(A,1) != size(A,2)
	error("factorMUMPS: Matrix must be square!")
    end
    n  = size(A,1);
    mumpsstat = [0];
    tic()
    p  = ccall( (:factor_mumps_cmplx_, MUMPSlibPath),
    		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
                 &n, &sym, &ooc,  convert(Ptr{Complex128}, pointer(A.nzval)), A.rowval, A.colptr, mumpsstat);
    checkMUMPSerror(mumpsstat);
    facTime = toq()
    factor = MUMPSfactorization(p,myid(),n,A.nzval[1],facTime);
    return factor
end

function factorMUMPS(A::SparseMatrixCSC{Float64},sym::Int=0,ooc::Int=0)
    # Generate LU-factorization of a real matrix A, 'sym' can be: 0=unsymmetric, 1=symm. pos def, 2=general symmetric
    if size(A,1) != size(A,2)
		error("factorMUMPS: Matrix must be square!")
    end
    n  = size(A,1);
    mumpsstat = [0];
    tic()
    p  = ccall( (:factor_mumps_, MUMPSlibPath),
 	       Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
               &n, &sym, &ooc, A.nzval, A.rowval, A.colptr, mumpsstat);
    checkMUMPSerror(mumpsstat);
    facTime = toq()
    factor = MUMPSfactorization(p,myid(),n,A.nzval[1],facTime);
    return factor
end


function checkMUMPSerror(mumpsstat)
	# check if MUMPS reported error
	
	if mumpsstat[1]==-10
	     error("MUMPS: Numerically singular matrix.");
	elseif mumpsstat[1]==-13
	     error("MUMPS: memory allocation error")
	elseif mumpsstat[1]==-40
	     error("MUMPS: matrix is not positive definite")
	elseif mumpsstat[1]==-90
	     error("MUMPS: Error in out-of-core management.Probably there is not enough disk space.")
	elseif mumpsstat[1]<0
	     error( @sprintf("MUMPS: error --> %d <--. Please refer to Ch. 7 of MUMPS User's guide!",mumpsstat[1]))
	end
end

function applyMUMPS{T1,T2,N}(factor::MUMPSfactorization{T1},rhs::AbstractArray{T2,N},x::Array=zeros(promote_type(T1,T2),size(rhs)),tr::Int=0)

	id1 = myid(); id2 = factor.worker
	if id1 != id2
		warn("Worker $id1 has no access to MUMPS factorization stored on $id2. Trying to remotecall!")
		return remotecall_fetch(factor.worker,applyMUMPS,factor,rhs,x,tr)::Array{promote_type(T1,T2),N}
	end
	 
	if size(rhs,1) != factor.n;  
		error("applyMUMPS: wrong size of rhs, size(A)=$(factor.n), size(rhs)=$(size(rhs,1)) x $(size(rhs,2))."); 
	end
    
	if size(x)!=size(rhs); 
		if isempty(x)
			x = zeros(promote_type(T1,T2),size(rhs))
		else
			error("applyMUMPS: wrong size of x, size(A)=$(factor.n), size(rhs)=$(size(rhs)), size(x)=$(size(x)) provided"); 
		end
	end
	
	return applyMUMPS!(factor,rhs,x,tr)::Array{promote_type(T1,T2),N}
end

function applyMUMPS!(factor::MUMPSfactorization{Float64},rhs::Array{Float64},
                      x::Array{Float64},tr::Int=0)

	 nrhs = size(rhs,2)
	ptr = factor.ptr
	ccall( (:solve_mumps_, MUMPSlibPath),
				Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64} ),
							&ptr,      &nrhs, 	 	rhs,      x, &tr)
	return x
end

function applyMUMPS!(factor::MUMPSfactorization{Float64},rhs::SparseMatrixCSC{Float64}, x::Array{Float64,2},tr::Int=0)
	nrhs = size(rhs,2)
	nzrhs = nnz(rhs)
	ptr = factor.ptr
	ccall( (:solve_mumps_sparse_rhs_, MUMPSlibPath),
              Void, (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64},
	      	  Ptr{Float64}, Ptr{Int64} ),
		      &ptr, &nzrhs, &nrhs, rhs.nzval, rhs.rowval, rhs.colptr, x, &tr)

	return x
end

function applyMUMPS!(factor::MUMPSfactorization{Complex128},rhs::Array{Complex128},
                      x::Array{Complex128},tr::Int=0)
	n     = size(rhs,1)
	nrhs  = size(rhs,2)
	ptr = factor.ptr
	ccall( (:solve_mumps_cmplx_, MUMPSlibPath),
				Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex64}, Ptr{Int64} ),
				&ptr,        &nrhs,   convert(Ptr{Complex128}, pointer(rhs)),   convert(Ptr{Complex128}, pointer(x)),   &tr)
	return x
end

function applyMUMPS!(factor::MUMPSfactorization{Complex128},rhs::SparseMatrixCSC{Complex128},
                      x::Array{Complex128,2},tr::Int=0)
nrhs = size(rhs,2)
nzrhs = nnz(rhs)
ptr = factor.ptr
ccall( (:solve_mumps_cmplx_sparse_rhs_, MUMPSlibPath),
			Void, (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Int64}, 
			Ptr{Int64}, Ptr{Complex64}, Ptr{Int64} ),
			&ptr, &nzrhs,&nrhs, convert(Ptr{Complex128}, pointer(rhs.nzval)), rhs.rowval, 
			rhs.colptr, convert(Ptr{Complex128}, pointer(x)), &tr)
return x
end

	

function destroyMUMPS(factor::MUMPSfactorization{Float64})
	 #  free memory
	id1 = myid(); id2 = factor.worker;
	if myid() != factor.worker
		remotecall_fetch(id2,destroyMUMPS,factor)
		return
	end
	ccall( (:destroy_mumps_, MUMPSlibPath),
			Int64, (Ptr{Int64}, ), &factor.ptr )
	factor.ptr = -1
	factor.n    = -1
	factor.time = -1.0

end

function destroyMUMPS(factor::MUMPSfactorization{Complex128})
	 #  free memory
	id1 = myid(); id2 = factor.worker;
	if myid() != factor.worker
		remotecall_fetch(id2,destroyMUMPS,factor)
		return
	end
	 ccall( (:destroy_mumps_cmplx_, MUMPSlibPath),
	 		Int64, (Ptr{Int64}, ), &factor.ptr )
	factor.ptr = -1
	factor.n    = -1
	factor.time = -1.0

end
