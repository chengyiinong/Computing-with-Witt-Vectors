import AbstractAlgebra

struct WittVectorsFq
    base_ring::AbstractAlgebra.ring
    precision::Int
    sum::Vector
    product::Vector
end

struct WittVectorsFqElement
    parent::WittvectorsFq
    elements::Vector
end

# To map W(F_p) to Z_p, for every Witt vectors (X_0, X_1, X_2, ...),
# we map to the p-adic integer w(X_0) + w(X_1) * p + w(X_2) * p ** 2 + ...
# where w(X_i) is implemented using the method teichmuller(X_i).
# For the inverse, for each p-adic integer a_0 + a_1 * p + a_2 * p ** 2 + ... ,
# reduce each a_i mod p, then the corresponding Witt vector is (a_0 mod p, a_1 mod p, ...)

function WittVectorsToZqElement(R::AbstractAlgebra.ring, n::Int, p::Int = Int(characteristic(R)))
    if 