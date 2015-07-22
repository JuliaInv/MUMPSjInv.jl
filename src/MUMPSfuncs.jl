function solveMUMPS(A::SparseMatrixCSC, rhs::arrayOrSparseCSC, x::Array=[], sym=0,ooc=0,tr=0)

	n    = size(rhs,1)
	nrhs = size(rhs,2)
	
	if isempty(x);
		x = zeros(eltype(A),n,nrhs)
	else
		if size(x)!=(n,nrhs); 
			error("applyMUMPS: wrong size of x provided")
		end
		# make x complex if necessary
		x = isreal(A) ? x : complex(x)
	end
	if (typeof(rhs) <: Array)
	  if norm(rhs)==0
		x = zeros(eltype(rhs),size(rhs))
		return x
	  end
	elseif (typeof(rhs) <: SparseMatrixCSC)
	  if nnz(rhs)==0
		x = zeros(eltype(rhs),size(rhs))
		return x
	  end
	end
	
	# factorization
	factor = factorMUMPS(A,sym,ooc)
	
	# solve system
	x = applyMUMPS(factor,rhs,x,tr)
	
	# free memory
	destroyMUMPS(factor)
        return x
end

function factorMUMPS(A::SparseMatrixCSC{Complex128,Int},sym=0,ooc=0)
    # Generate LU-factorization of complex matrix A, 'sym' can be: 0=unsymmetric, 1=symm. pos def, 2=general symmetric
    if size(A,1) != size(A,2)
	error("factorMUMPS: Matrix must be square!")
    end
    n  = size(A,1);
    mumpsstat = [0];
    tic()
    p  = ccall( (:factor_mumps_cmplx_, "/home/patrick/Julia/v0.4/MUMPS.jl/lib/MUMPS"),
    		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
                 &n, &sym, &ooc,  convert(Ptr{Complex128}, pointer(A.nzval)), A.rowval, A.colptr, mumpsstat);
    checkMUMPSerror(mumpsstat);
    facTime = toq()
    factor = MUMPSfactorizationComplex(p,myid(),n,facTime);
    return factor
end

function factorMUMPS(A::SparseMatrixCSC{Float64,Int},sym=0,ooc=0)
    # Generate LU-factorization of a real matrix A, 'sym' can be: 0=unsymmetric, 1=symm. pos def, 2=general symmetric
    if size(A,1) != size(A,2)
		error("factorMUMPS: Matrix must be square!")
    end
    n  = size(A,1);
    mumpsstat = [0];
    tic()
    p  = ccall( (:factor_mumps_, "/home/patrick/Julia/v0.4/MUMPS.jl/lib/MUMPS"),
 	       Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
               &n, &sym, &ooc, A.nzval, A.rowval, A.colptr, mumpsstat);
    checkMUMPSerror(mumpsstat);
    facTime = toq()
    factor = MUMPSfactorizationReal(p,myid(),n,facTime);
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

function applyMUMPS(factor::MUMPSfactorizationReal,rhs::arrayOrSparseCSCReal,
                      x::Array{Float64}=zeros(0,0),tr=0)
	# solve for right hand side(s)

	id1 = myid(); id2 = factor.worker
	if id1 != id2
		warn("Worker $id1 has no access to MUMPS factorization stored on $id2. Trying to remotecall!")
		return remotecall_fetch(factor.worker,applyMUMPS,factor,rhs,x,tr)::Array{Float64,2}
	end
	
	 n    = size(rhs,1)
	 nrhs = size(rhs,2)
	
	# check size of rhs, allocate space for x
	if n != factor.n;  error("applyMUMPS: wrong size of rhs, size(A)=$(factor.n), size(rhs)=$n x $nrhs."); end
	x = isempty(x) ? x = zeros(n,nrhs) : x
	if size(x)!=(n,nrhs); error("applyMUMPS: wrong size of x provided"); end
	
	ptr = factor.ptr
	if typeof(rhs) <: Array
		ccall( (:solve_mumps_, "/home/patrick/Julia/v0.4/MUMPS.jl/lib/MUMPS"),
				Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64} ),
							&ptr,      &nrhs, 	 	rhs, x, &tr)
	else
		nzrhs = nnz(rhs)
		ccall( (:solve_mumps_sparse_rhs_, "/home/patrick/Julia/v0.4/MUMPS.jl/lib/MUMPS"),
	              Void, (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64},
		      Ptr{Float64}, Ptr{Int64} ),
		      &ptr, &nzrhs, &nrhs, rhs.nzval, rhs.rowval, rhs.colptr, x, &tr)
	end
	return x::Array{Float64,2}
end

function applyMUMPS(factor::MUMPSfactorizationComplex,rhs::arrayOrSparseCSCComplex,
                      x::Array{Complex{Float64}}=zeros(Complex{Float64},0,0),tr=0)
	# solve for right hand side(s)

	id1 = myid(); id2 = factor.worker
	if id1 != id2
		warn("Worker $id1 has no access to MUMPS factorization stored on $id2. Trying to remotecall!")
		return remotecall_fetch(factor.worker,applyMUMPS,factor,rhs,x,tr)::Array{Complex{Float64},2}
	end
	
	n     = size(rhs,1)
	nrhs  = size(rhs,2)
	
	# check size of rhs, allocate space for x
	if n != factor.n;  error("applyMUMPS: wrong size of rhs, size(A)=$(factor.n), size(rhs)=$n x $nrhs."); end
	x = isempty(x) ? x = zeros(Complex{Float64},n,nrhs) : x
	if size(x)!=(n,nrhs); error("applyMUMPS: wrong size of x provided"); end
	
	ptr = factor.ptr
	if typeof(rhs) <: Array
		ccall( (:solve_mumps_cmplx_, "/home/patrick/Julia/v0.4/MUMPS.jl/lib/MUMPS"),
					Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex64}, Ptr{Int64} ),
					&ptr,        &nrhs,   convert(Ptr{Complex128}, pointer(rhs)),   convert(Ptr{Complex128}, pointer(x)),   &tr)
	else
		nzrhs = nnz(rhs)
		ccall( (:solve_mumps_cmplx_sparse_rhs_, "/home/patrick/Julia/v0.4/MUMPS.jl/lib/MUMPS"),
					Void, (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Int64}, 
					Ptr{Int64}, Ptr{Complex64}, Ptr{Int64} ),
					&ptr, &nzrhs,&nrhs, convert(Ptr{Complex128}, pointer(rhs.nzval)), rhs.rowval, 
					rhs.colptr, convert(Ptr{Complex128}, pointer(x)), &tr)
	end
	return x::Array{Complex{Float64},2}
end

function destroyMUMPS(factor::MUMPSfactorizationReal)
	 #  free memory
	id1 = myid(); id2 = factor.worker;
	if myid() != factor.worker
		remotecall_fetch(id2,destroyMUMPS,factor)
		return
	end
	ccall( (:destroy_mumps_, "/home/patrick/Julia/v0.4/MUMPS.jl/lib/MUMPS"),
			Int64, (Ptr{Int64}, ), &factor.ptr )
	factor.ptr = -1
	factor.n    = -1
	factor.time = -1.0

end

function destroyMUMPS(factor::MUMPSfactorizationComplex)
	 #  free memory
	id1 = myid(); id2 = factor.worker;
	if myid() != factor.worker
		remotecall_fetch(id2,destroyMUMPS,factor)
		return
	end
	 ccall( (:destroy_mumps_cmplx_, "/home/patrick/Julia/v0.4/MUMPS.jl/lib/MUMPS"),
	 		Int64, (Ptr{Int64}, ), &factor.ptr )
	factor.ptr = -1
	factor.n    = -1
	factor.time = -1.0

end
