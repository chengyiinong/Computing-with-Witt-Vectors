# Computing with Witt Vectors

This repository contains a Julia code implementation `Twoples_WittVectors.jl` on the topic of computing Witt vectors $`W(\mathbb{F}_{p})`$ and $`W(\mathbb{F}_{p^f})`$, which is part of a final project for a Twoples Reading Project (2025 Fall), titled "Non-Archimedean Banach Algebras and Adic Spaces", mentored by Gabriel Ong. Using the computer algebra package `Oscar`, the Witt vectors $`W(\mathbb{F}_{p})`$ and $`W(\mathbb{F}_{p^f})`$ are defined, and the operations of addition, subtraction, multiplication, Fronebii, and Verschiebungen for Witt vectors are implemented in the Julia language. 

Despite being able to implement $`W(\mathbb{F}_{p})`$ and the operations on $`W(\mathbb{F}_{p})`$ fully, there are limitations in the current version of `Twoples_WittVectors.jl` in defining the Witt vectors $`W(\mathbb{F}_{p^f})`$ for elements which are not in the base field $`\mathbb{F}_p`$. Further changes also have to be made to implement the isomorphism between $`W(\mathbb{F}_{p^f})`$ and $`\mathbb{Z}_p[\zeta_{p^f-1}]`$. 

At last, a sample test case is implemented in `test_cases.jl` for $W(\mathbb{F}_{13})$. 
