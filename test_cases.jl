# Note: The code doesn't yet support operations on elements from F_{p^f} which are not in the base field F_p  

include("Twoples_WittVectors.jl")

F, a = finite_field(13, 1, "a")

witt_ring = WittVectorsFq(F, 100)

U = WittVectorsFqElement{FinFieldElem}(witt_ring, [F(7), F(4), F(2), F(5), F(1), F(8), F(12), F(10)])

W = WittVectorsFqElement{FinFieldElem}(witt_ring, [F(3), F(4), F(2), F(5), F(7), F(11), F(2)])

wittvector_U = WittVector(witt_ring, U)

wittvector_W = WittVector(witt_ring, W)

sum = plus(witt_ring, U, W)

sub = subtract(witt_ring, U, W)

product = multiply(witt_ring, U, W)

Frobenii_U = Frobenii(witt_ring, U)

Frobenii_W = Frobenii(witt_ring, W)

Verschiebungen_U = Verschiebungen(witt_ring, U)

Verschiebungen_W = Verschiebungen(witt_ring, W)

println("Witt vector U = ", wittvector_U, "\n")

println("Witt vector W = ", wittvector_W, "\n")

println("U + W = ", sum, "\n")

println("U - W = ", sub, "\n")

println("U * W = ", product, "\n")

println("Frobenii(U)", Frobenii_U, "\n")

println("Frobenii(W)", Frobenii_W, "\n")

println("Verschiebungen(U)", Verschiebungen_U, "\n")

println("Verschiebungen(W)", Verschiebungen_W, "\n")
