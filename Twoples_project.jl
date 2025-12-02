import AbstractAlgebra
using Oscar

# This defines the datatype for the ring of Witt vectors, where users can input the 
# base ring and the precision used

struct WittVectorsFq
    base_ring::Oscar.Ring
    qpower::Int
    precision::Int
end

# This defines the elements of a Witt vector in a vector type

struct WittVectorsFqElement{T}
    parent::WittVectorsFq
    elements::Vector{T}
end

# This function fills up the unspecified entries of a Witt vector with 0 until the precision length 

function WittVector(F::WittVectorsFq, W::WittVectorsFqElement)
    if length(W.elements) < R.precision 
        zero_arr = fill(ZZRingElem(0), R.precision - length(W.elements))
        new_element = vcat(W.elements, zero_arr)
    else
        new_element = W.elements
    end
    return new_element
end

# To map W(F_p) to Z_p, for every Witt vectors (X_0, X_1, X_2, ...),
# we map to the p-adic integer w(X_0) + w(X_1) * p + w(X_2) * p ** 2 + ...
# where w(X_i) is implemented using the method teichmuller(X_i).
# For the inverse, for each p-adic integer a_0 + a_1 * p + a_2 * p ** 2 + ... ,
# reduce each a_i mod p, then the corresponding Witt vector is (a_0 mod p, a_1 mod p, ...)

function WittVectorsToZqElement(F::WittVectorsFq, W::WittVectorsFqElement)
    R = F.base_ring
    prec = F.precision
    p = prime(R)
    padic_eq = R(0)
    f = F.qpower
    vec = WittVector(F, W)
    if f == 1
        for i in 1:prec
            qadic_eq += R(teichmuller(R(vec[i])) * p^(i - 1))
        end
    else 
        for i in 1:prec
            r = mod(i, f)
            qadic_eq += R((teichmuller(R(vec[i])) ** (p ** (f - r))) * p^(i - 1))
        end
    return qadic_eq
end
