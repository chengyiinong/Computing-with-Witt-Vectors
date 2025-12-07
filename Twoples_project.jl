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

function WittVector(R::WittVectorsFq, W::WittVectorsFqElement)
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
            padic_eq += R(teichmuller(R(vec[i])) * p^(i - 1))
        end
    return padic_eq
    else 
        for i in 1:prec
            r = mod(i, f)
            padic_eq += R((teichmuller(R(vec[i]))^(p^(f - r))) * p^(i - 1))
        end
    return padic_eq
    end
end

# The problem to solve: It's very hard to isolate out the leading coefficient of a p-adic number
# The only way I have now is to just run a for loop by keep on subtracting 1 until you hit 0 for the leading coeff.
# !!!The solution is to use p = lift(Zx, f) where Zx, x = ZZ["x"] and f is the p-adic number   
# For unramified extensions, consider looking into the cyclotomic integers and find ways to 
# realize the teichmuller map from W(Fq) to Zp[zeta_{q-1}]

function plus(F::WittVectorsFq, X::WittVectorsFqElement, Y::WittVectorsFqElement)
    R = F.base_ring
    prec = F.precision
    p = prime(R)
    f = F.qpower
    Zx, x = ZZ["x"]
    qadic_X = WittVectorsToZqElement(F, X)
    qadic_Y = WittVectorsToZqElement(F, Y)
    sum = R(qadic_X + qadic_Y)
    Z0 = (X[1] + Y[1]) % p^f
    Z_elements = []
    push!(Z_elements, Z0)
    if f == 1
        for i in 1:prec
            num = R(0)
            for j in 1:i 
                num += R(teichmuller(R(Z_elements[j])) * p^(j - 1))
            end
            rem = sum - num
            Zi = lift(Zx, rem)(1) % p^i
            push!(Z_elements, Zi)
        end
        return Z_elements
    else
        for i in 1:prec
            num = R(0)
            for j in 1:i 
                r = mod(j, f)
                num += R(teichmuller(R(Z_elements[j])^(p^(f - r))) * p^(j - 1))
            end
            rem = sum - num
            Zi = ((lift(Zx, rem)(1) % p^i)^p^i) % p^f
            push!(Z_elements, Zi)
        end
        return Z_elements
    end
end

function subtract(F::WittVectorsFq, X::WittVectorsFqElement, Y::WittVectorsFqElement)
    R = F.base_ring
    prec = F.precision
    p = prime(R)
    f = F.qpower
    Zx, x = ZZ["x"]
    qadic_X = WittVectorsToZqElement(F, X)
    qadic_Y = WittVectorsToZqElement(F, Y)
    sum = R(qadic_X - qadic_Y)
    Z0 = (X[1] - Y[1]) % p^f
    Z_elements = []
    push!(Z_elements, Z0)
    if f == 1
        for i in 1:prec
            num = R(0)
            for j in 1:i 
                num += R(teichmuller(R(Z_elements[j])) * p^(j - 1))
            end
            rem = sum - num
            Zi = lift(Zx, rem)(1) % p^i
            push!(Z_elements, Zi)
        end
        return Z_elements
    else
        for i in 1:prec
            num = R(0)
            for j in 1:i 
                r = mod(j, f)
                num += R(teichmuller(R(Z_elements[j])^(p^(f - r))) * p^(j - 1))
            end
            rem = sum - num
            Zi = ((lift(Zx, rem)(1) % p^i)^(p^i)) % p^f
            push!(Z_elements, Zi)
        end
        return Z_elements
    end
end

function multiply(F::WittVectorsFq, X::WittVectorsFqElement, Y::WittVectorsFqElement)
    R = F.base_ring
    prec = F.precision
    p = prime(R)
    f = F.qpower
    Zx, x = ZZ["x"]
    qadic_X = WittVectorsToZqElement(F, X)
    qadic_Y = WittVectorsToZqElement(F, Y)
    prod = R(qadic_X * qadic_Y)
    Z0 = (X[1] * Y[1]) % p^f
    Z_elements = []
    push!(Z_elements, Z0)
    if f == 1
        for i in 1:prec
            num = R(0)
            for j in 1:i 
                num += R(teichmuller(R(Z_elements[j])) * p^(j - 1))
            end
            rem = prod - num
            Zi = lift(Zx, rem)(1) % p^i
            push!(Z_elements, Zi)
        end
        return Z_elements
    else
        for i in 1:prec
            num = R(0)
            for j in 1:i 
                r = mod(j, f)
                num += R(teichmuller(R(Z_elements[j])^(p^(f - r))) * p^(j - 1))
            end
            rem = prod - num
            Zi = ((lift(Zx, rem)(1) % p^i)^(p^i)) % p^f
            push!(Z_elements, Zi)
        end
        return Z_elements
    end
end
            
function Frobenii(F::WittVectorsFq, W::WittVectorsFqElement)
    vec = WittVector(F, W)
    prec = F.precision
    p = prime(F.base_ring)
    for i in 1:prec
        vec = vec^p
    end
    return vec 
end

function Verschiebungen(F::WittVectorsFq, W::WittVectorsFqElement)
    vec = WittVector(F, W)
    zero_arr = [0]
    img_vec = vcat(zero_arr, vec)
    return img_vec
end
