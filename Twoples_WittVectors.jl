import AbstractAlgebra
using Oscar

# This defines the datatype for the ring of Witt vectors, where users can input the 
# base ring (i.e., $F_p$ or $F_{p^f}$) and the precision used

struct WittVectorsFq
    base_ring::Oscar.FinField
    precision::Int
end

# This defines the elements of a Witt vector in a vector type

struct WittVectorsFqElement{T}
    parent::WittVectorsFq
    elements::Vector{T}
end

@doc raw"""
    WittVector(R::WittVectorsFq, W::WittVectorsFqElement) -> vector{FinFieldElem}

Given a Witt vector of $l <= n$ specified elements $F_p$ or $F_{p^k}$, return a Witt
vector of precision $n$ by filling up the $n - l$ unspecified elements with the zero element 0.
"""
function WittVector(R::WittVectorsFq, W::WittVectorsFqElement)
    if length(W.elements) < R.precision 
        zero_arr = fill(R.base_ring(0), R.precision - length(W.elements))
        new_element = vcat(W.elements, zero_arr)
    else
        new_element = W.elements
    end
    return new_element
end

@doc raw"""
    WittVectorsToZqElement(F::WittVectorsFq, W::WittVectorsFqElement) -> QadicFieldElem

Return the corresponding $q$-adic representation for a given Witt vector $(X_1, X_2, ..., X_n, ...)$
under the isomorphism between Witt vectors ($W(F_p)$ or $W(F_{p^f})$) and q-adic extensions ($Z_p$ or $Z_p[zeta_{p^f - 1}]$). 
"""
function WittVectorsToZqElement(F::WittVectorsFq, W::WittVectorsFqElement)
    R = F.base_ring
    prec = F.precision
    f = degree(R)
    p = characteristic(prime_field(R))
    Q, q = qadic_field(p, f, precision = prec)
    padic_eq = Q(0)
    vec = WittVector(F, W)
    if f == 1
        for i in 1:prec
            padic_eq += teichmuller(Q(lift(ZZ, vec[i]))) * p^(i - 1)
        end
    return padic_eq
    else 
        for i in 1:prec
            r = mod(i, f)
            padic_eq += (teichmuller(Q(lift(ZZ, vec[i])))^(p^(f - r))) * p^(i - 1)
        end
    return padic_eq
    end
end

###############################################################################
#
#   Operations on the ring of Witt vectors  
#
###############################################################################

@doc raw"""
    plus(F::WittVectorsFq, X::WittVectorsFqElement, Y::WittVectorsFqElement) -> vector{FinFieldElem}

Return the sum of two Witt vectors $X$ and $Y$ as a vector of $Z_p$ or $Z_p[zeta_{p^f - 1}]$
through first mapping $X$ and $Y$ to their corresponding $q$-adic representation, then adding $Z = X + Y$ in the 
$q$-adic ring, and mapping $Z$ back to its corresponding Witt vector. 
"""
function plus(F::WittVectorsFq, X::WittVectorsFqElement, Y::WittVectorsFqElement)
    R = F.base_ring
    prec = F.precision
    f = degree(R)
    p = characteristic(prime_field(R))
    Zx, x = ZZ["x"]
    Q, q = qadic_field(p, f, precision = prec)
    qadic_X = WittVectorsToZqElement(F, X)
    qadic_Y = WittVectorsToZqElement(F, Y)
    sum = Q(qadic_X + qadic_Y)
    Z0 = WittVector(F, X)[1] + WittVector(F, Y)[1]
    Z_elements = []
    push!(Z_elements, Z0)
    if f == 1
        for i in 1:prec
            num = Q(0)
            for j in 1:i 
                num += teichmuller(Q(lift(ZZ, Z_elements[j]))) * p^(j - 1)
            end
            rem = (sum - num) / p^i
            Zi = R(lift(Zx, rem)(1) % p)
            push!(Z_elements, Zi)
        end
        return Z_elements
    else
        for i in 1:prec
            num = Q(0)
            for j in 1:i 
                r = mod(j, f)
                num += teichmuller(Q(lift(ZZ, Z_elements[j]))^(p^(f - r))) * p^(j - 1)
            end
            rem = (sum - num) / p^i
            Zi = R((lift(Zx, rem)(1) % p)^p^i)
            push!(Z_elements, Zi)
        end
        return Z_elements
    end
end

@doc raw"""
    subtract(F::WittVectorsFq, X::WittVectorsFqElement, Y::WittVectorsFqElement) -> vector{FinFieldElem}

Return the subtraction of two Witt vectors $X$ and $Y$ as a vector of $Z_p$ or $Z_p[zeta_{p^f - 1}]$
through first mapping $X$ and $Y$ to their corresponding $q$-adic representation, then subtracting $Z = X - Y$ in the 
$q$-adic ring, and mapping $Z$ back to its corresponding Witt vector. 
"""
function subtract(F::WittVectorsFq, X::WittVectorsFqElement, Y::WittVectorsFqElement)
    R = F.base_ring
    prec = F.precision
    f = degree(R)
    p = characteristic(prime_field(R))
    Zx, x = ZZ["x"]
    Q, q = qadic_field(p, f, precision = prec)
    qadic_X = WittVectorsToZqElement(F, X)
    qadic_Y = WittVectorsToZqElement(F, Y)
    sum = Q(qadic_X - qadic_Y)
    Z0 = WittVector(F, X)[1] - WittVector(F, Y)[1]
    Z_elements = []
    push!(Z_elements, Z0)
    if f == 1
        for i in 1:prec
            num = Q(0)
            for j in 1:i 
                num += teichmuller(Q(lift(ZZ, Z_elements[j]))) * p^(j - 1)
            end
            rem = (sum - num) / p^i
            Zi = R(lift(Zx, rem)(1) % p)
            push!(Z_elements, Zi)
        end
        return Z_elements
    else
        for i in 1:prec
            num = Q(0)
            for j in 1:i 
                r = mod(j, f)
                num += teichmuller(Q(lift(ZZ, Z_elements[j]))^(p^(f - r))) * p^(j - 1)
            end
            rem = (sum - num) / p^i
            Zi = R((lift(Zx, rem)(1) % p)^p^i)
            push!(Z_elements, Zi)
        end
        return Z_elements
    end
end

@doc raw"""
    multiply(F::WittVectorsFq, X::WittVectorsFqElement, Y::WittVectorsFqElement) -> vector{FinFieldElem}

Return the multiplication of two Witt vectors $X$ and $Y$ as a vector of $Z_p$ or $Z_p[zeta_{p^f - 1}]$
through first mapping $X$ and $Y$ to their corresponding $q$-adic representation, then multiplying $Z = X * Y$ in the 
$q$-adic ring, and mapping $Z$ back to its corresponding Witt vector. 
"""
function multiply(F::WittVectorsFq, X::WittVectorsFqElement, Y::WittVectorsFqElement)
    R = F.base_ring
    prec = F.precision
    f = degree(R)
    p = characteristic(prime_field(R))
    Zx, x = ZZ["x"]
    Q, q = qadic_field(p, f, precision = prec)
    qadic_X = WittVectorsToZqElement(F, X)
    qadic_Y = WittVectorsToZqElement(F, Y)
    sum = Q(qadic_X * qadic_Y)
    Z0 = WittVector(F, X)[1] * WittVector(F, Y)[1]
    Z_elements = []
    push!(Z_elements, Z0)
    if f == 1
        for i in 1:prec
            num = Q(0)
            for j in 1:i 
                num += teichmuller(Q(lift(ZZ, Z_elements[j]))) * p^(j - 1)
            end
            rem = (sum - num) / p^i
            Zi = R(lift(Zx, rem)(1) % p)
            push!(Z_elements, Zi)
        end
        return Z_elements
    else
        for i in 1:prec
            num = Q(0)
            for j in 1:i 
                r = mod(j, f)
                num += teichmuller(Q(lift(ZZ, Z_elements[j]))^(p^(f - r))) * p^(j - 1)
            end
            rem = (sum - num) / p^i
            Zi = R((lift(Zx, rem)(1) % p)^p^i)
            push!(Z_elements, Zi)
        end
        return Z_elements
    end
end

@doc raw"""
    Frobenii(F::WittVectorsFq, W::WittVectorsFqElement) -> vector{FinFieldElem}

Given a Witt vector $W = (X_0, X_1, X_2, ...)$, return the image of $W$ under the 
Frobenii map, i.e., $(X_0^p, X_1^p, X_2^p, ...)$.
"""
function Frobenii(F::WittVectorsFq, W::WittVectorsFqElement)
    R = F.base_ring
    vec = WittVector(F, W)
    prec = F.precision
    p = characteristic(prime_field(R))
    for i in 1:prec
        vec = vec^p
    end
    return vec 
end

@doc raw"""
    Verschiebungen(F::WittVectorsFq, W::WittVectorsFqElement) -> vector{FinFieldElem}

Given a Witt vector $W = (X_0, X_1, X_2, ...)$, return the image of $W$ under the 
Verschiebungen map, i.e., $(0, X_0, X_1, X_2, ...)$.
"""
function Verschiebungen(F::WittVectorsFq, W::WittVectorsFqElement)
    R = F.base_ring
    vec = WittVector(F, W)
    zero_arr = [R(0)]
    img_vec = vcat(zero_arr, vec)
    return img_vec
end
