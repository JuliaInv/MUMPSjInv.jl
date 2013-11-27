function solveMUMPS(A::SparseMatrixCSC, rhs::Array, x::Array=[], tr=0, sym=1)

	n    = size(rhs,1)
	nrhs = size(rhs,2)
	ncol = size(A,2)
	nrow = size(A,1)
	
	# check size of rhs, allocate space for x
	if tr==1
		if n != ncol;  error("solveMUMPS: wrong size of rhs"); end
		x = isempty(x) ? x = zeros(nrow,nrhs) : x
		if size(x)!=(nrow,nrhs); error("applyMUMPS: wrong size of x provided");end
    else
		if n != nrow; error("solveMUMPS: wrong size of rhs"); end
		x = isempty(x) ? x = zeros(ncol,nrhs) : x
		if size(x)!=(ncol,nrhs); error("applyMUMPS: wrong size of x provided");end
    end

	# make x complex if necessary
    x = isreal(rhs) ? x : complex(x)

	# factorization
	factor = factorMUMPS(A,sym)
	
	# solve system
	x = applyMUMPS(factor,rhs,x,tr)
	
	# free memory
	destroyMUMPS(factor)

   return x
end


function factorMUMPS(A::SparseMatrixCSC,sym=1)
	# Generate LU-factorization of real matrix A, 'sym' can be: 0=unsymmetric, 1=symm. pos def, 2=general symmetric
	
    nrow  = size(A,1);
    ncol  = size(A,2);
    mumpsstat = [0];
	isReal = iseltype(A,Real)
    
    tic();
		if isReal
    		p  = ccall( (:factor_mumps_, "../lib/MUMPS"),
    			 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
    	         &nrow, &sym, A.nzval, A.rowval, A.colptr, mumpsstat);
    	else
		    p  = ccall( (:factor_mumps_cmplx_, "../lib/MUMPS"),
		    		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Complex64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
		             &nrow, &sym, A.nzval, A.rowval, A.colptr, mumpsstat);
		end
    facTime = toc();
    
    checkMUMPSerror(mumpsstat);
    
    factor = MUMPSfactorization(p,nrow,ncol,isReal,facTime);
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
	elseif mumpsstat[1]<0
	     error( @sprintf("MUMPS: error --> %d <--. Please refer to Ch. 7 of MUMPS User's guide!",mumpsstat[1]))
	end
end

function applyMUMPS(factor::MUMPSfactorization,rhs::Array,x::Array=[],tr=0)
	# solve for right hand side(s)
	
    n    = size(rhs,1)
    nrhs = size(rhs,2)
    rhs  = reshape(rhs,n,nrhs)

	# check size of rhs, allocate space for x
	if tr==1
		if n != factor.ncol;  error("applyMUMPS: wrong size of rhs"); end
		x = isempty(x) ? x = zeros(factor.nrow,nrhs) : x
        if size(x)!=(factor.nrow,nrhs); error("applyMUMPS: wrong size of x provided"); end
    else
		if n != factor.nrow; error("applyMUMPS: wrong size of rhs"); end
		x = isempty(x) ? x = zeros(factor.ncol,nrhs) : x
        if size(x)!=(factor.ncol,nrhs); error("applyMUMPS: wrong size of x provided"); end
    end
    
	if isreal(rhs)!=factor.real
	   factor.real  ? error("applyMUMPS: rhs must be real.") : 	error("applyMUMPS: rhs must be complex.")
	end

	 # make x complex if necessary
	 x = factor.real ? x : complex(x)
    
     ptr = factor.ptr
	 if factor.real
     	ccall( (:solve_mumps_, "../lib/MUMPS"),
            Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64} ),
                        &ptr,      &nrhs, 	 	rhs, x, &tr)
	else
		ccall( (:solve_mumps_cmplx_, "../lib/MUMPS"),
		       Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Complex64}, Ptr{Complex64}, Ptr{Int64} ),
		               &ptr,        &nrhs,     rhs,            x,              &tr)
	end
    return x
end

function destroyMUMPS(factor::MUMPSfactorization)
    #  free memory
        ptr = factor.ptr
    ccall( (:destroy_mumps_, "../lib/MUMPS"),
               Int64, (Ptr{Int64}, ), &ptr );
end
