function solveMUMPS(A::SparseMatrixCSC, rhs::Array, x::Array=[], sym=0,tr=0)

	n    = size(rhs,1)
	nrhs = size(rhs,2)
	
	if size(A,1)!=n || size(A,2) != n;
		error("solveMUMPS: matrix must be square and match length(rhs).")
	end
	
	if isempty(x);
		x = zeros(n,nrhs)
	elseif size(x)!=(n,nrhs); 
		error("applyMUMPS: wrong size of x provided")
	end
	
	if norm(rhs)==0
		x = rhs
		return x
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


function factorMUMPS(A::SparseMatrixCSC,sym=0)
	# Generate LU-factorization of real matrix A, 'sym' can be: 0=unsymmetric, 1=symm. pos def, 2=general symmetric
	
    if size(A,1) != size(A,2)
		error("factorMUMPS: Matrix must be square!")
  	end
	n  = size(A,1);
   	
 	mumpsstat = [0];
	isReal = iseltype(A,Real)
    
    tic();
		if isReal
    		p  = ccall( (:factor_mumps_, "../lib/MUMPS"),
    			 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
    	         &n, &sym, A.nzval, A.rowval, A.colptr, mumpsstat);
    	else
		    p  = ccall( (:factor_mumps_cmplx_, "../lib/MUMPS"),
		    		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Complex64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
		             &n, &sym,  convert(Ptr{Complex64}, pointer(A.nzval)), A.rowval, A.colptr, mumpsstat);
		end
    facTime = toc();
    
    checkMUMPSerror(mumpsstat);
    
    factor = MUMPSfactorization(p,n,isReal,facTime);
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
	if n != factor.n;  error("applyMUMPS: wrong size of rhs"); end
	x = isempty(x) ? x = zeros(n,nrhs) : x
    if size(x)!=(n,nrhs); error("applyMUMPS: wrong size of x provided"); end
    
	if factor.real && !isreal(rhs)
	   error("applyMUMPS: rhs must be real.")	
	end
	

	 # make x complex if necessary
	 x   = factor.real ? x   : complex(x)
     rhs = factor.real ? rhs : complex(rhs)
	
     ptr = factor.ptr
	 if factor.real
     	ccall( (:solve_mumps_, "../lib/MUMPS"),
            Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64} ),
                        &ptr,      &nrhs, 	 	rhs, x, &tr)
	else
		ccall( (:solve_mumps_cmplx_, "../lib/MUMPS"),
		       Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Complex64}, Ptr{Complex64}, Ptr{Int64} ),
		               &ptr,        &nrhs,   convert(Ptr{Complex64}, pointer(rhs)),   convert(Ptr{Complex64}, pointer(x)),   &tr)
	end
    return x
end

function destroyMUMPS(factor::MUMPSfactorization)
    #  free memory
	if factor.real
    	ccall( (:destroy_mumps_, "../lib/MUMPS"),
               Int64, (Ptr{Int64}, ), &factor.ptr );
	else
    	ccall( (:destroy_mumps_cmplx_, "../lib/MUMPS"),
               Int64, (Ptr{Int64}, ), &factor.ptr );
	end
end
