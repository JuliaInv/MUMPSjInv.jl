using MUMPS
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
include("getDivGrad.jl");
Ns = (8,16,24,32,48,64,72)
MUMPStime=zeros(length(Ns))
MUMPStimeOOC=zeros(length(Ns))
Juliatime=zeros(length(Ns))


for i=1:length(Ns)
	println(@sprintf("Solve div-grad on %d x %d x %d grid",Ns[i],Ns[i],Ns[i]));
	
	A = getDivGrad(Ns[i],Ns[i],Ns[i]);
	n = size(A,1);
    rhs = randn(n);
	
	# solve using mumps
	tic()
	x1 = solveMUMPS(A,rhs,1,1);
	MUMPStimeOOC[i]= toc();
	
	# solve using mumps in core
	tic()
	x2 = solveMUMPS(A,rhs,1,0);
	MUMPStime[i]= toc();

	# solve using mumps
	tic()
	x3 = A\rhs;
	Juliatime[i]= toc();
	
	@test norm(x1-x2)/norm(x1) < 1e-9
	@test norm(x1-x3)/norm(x1) < 1e-9
end
println(@sprintf("| %4s | %8s | %8s | %8s | %8s |","n","MUMPSooc","MUMPS","Julia","speedup"))

for i=1:length(Ns)
	println(@sprintf("| %d^3 | %1.5f | %1.5f | %1.5f | %1.5f |",Ns[i],MUMPStimeOOC[i],MUMPStime[i],Juliatime[i],Juliatime[i]/MUMPStime[i]))
end