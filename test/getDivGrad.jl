using SparseArrays
using LinearAlgebra

#----------------- Get the A matrix
function getDivGrad(n1,n2,n3)

        # the Divergence
        D1 = kron(sparse(1.0I, n3, n3),kron(sparse(1.0I, n2, n2),ddx(n1)))
        D2 = kron(sparse(1.0I, n3, n3),kron(ddx(n2),sparse(1.0I, n1, n1)))
        D3 = kron(ddx(n3),kron(sparse(1.0I, n2, n2),sparse(1.0I, n1, n1)))
        # DIV from faces to cell-centers
        Div = [D1 D2 D3]

        return Div*Div';
end
#----------------- 1D finite difference on staggered grid
function ddx(n)
# generate 1D derivatives
    return d = spdiagm(0 => -ones(n), 1 => ones(n))[1:end-1,:]
end
