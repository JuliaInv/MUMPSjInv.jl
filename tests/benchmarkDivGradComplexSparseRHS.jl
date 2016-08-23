using MUMPS
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
include("getDivGrad.jl");
Ns   = (32,48,64)
nrhs = (1,10,100,500)
rhsDensity = 0.01
MUMPStimeSparse =zeros(length(nrhs))
MUMPStimeDense  =zeros(length(nrhs))

for i=1:length(Ns)
	println(@sprintf("Factoring div-grad matrix on %d x %d x %d grid",Ns[i],Ns[i],Ns[i]))
	A = getDivGrad(Ns[i],Ns[i],Ns[i]);
	n = size(A,1);
	r = ones(n);
	A = A + im*spdiagm(r,0);
	
	#Factor the matrix
	Afac1 = factorMUMPS(A,1)
	Afac2 = factorMUMPS(A,1)
	
	rhsInit = spdiagm((ones(n),ones(n-Ns[i]),ones(n-2*Ns[i])),[0,Ns[i],2*Ns[i]],n,n)
	#rhsInit = sprandn(n,nrhs[end],rhsDensity)
	rhsInit = complex(rhsInit)
	for j = 1:length(nrhs)
	  println(@sprintf("Solve with %d rhs",nrhs[j]));
	  #rhs = sprandn(n,nrhs[j],rhsDensity)
	  rhs = rhsInit[:,1:nrhs[j]]
	  
	  # solve using mumps with sparse RHS
	  tic()
	  x1 = applyMUMPS(Afac1,rhs)
	  MUMPStimeSparse[j]= toc();
	
	  # solve using mumps with dense rhs
	  rhs = full(rhs)
	  tic()
	  x2 = applyMUMPS(Afac2,rhs)
	  MUMPStimeDense[j]= toc();
	  for k = 1:nrhs[j]
	    @test norm(x1[:,k]-x2[:,k])/norm(x1[:,k]) < 1e-9
	  end
        end
	destroyMUMPS(Afac1)
	destroyMUMPS(Afac2)

        println(@sprintf("| %4s | %8s | %8s | %8s |","nrhs","sparseRHS","denseRHS","speedup"))

        for j=1:length(nrhs)
	  println(@sprintf("| %d | %1.5f | %1.5f | %1.5f |",nrhs[j],MUMPStimeSparse[j],MUMPStimeDense[j],MUMPStimeDense[j]/MUMPStimeSparse[j]))
        end
        println(" ")
end
