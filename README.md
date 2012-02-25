hTensor
=======

A Haskell package for multidimensional arrays, simple tensor computations and multilinear algebra.

Array dimensions have an "identity" which is preserved in data manipulation. Indices are explicitly selected by name in expressions, and Einstein's summation convention for repeated indices is automatically applied.

The library has a purely functional interface: arrays are immutable, and operations work on whole structures which can be assembled and decomposed using simple primitives. Arguments are automatically made conformable by replicating them along extra dimensions appearing in an operation.

There is preliminary support for geometric algebra, multidimensional linear systems of equations, and tensor decompositions.

- [Source code and documentation][source]

- [Tutorial][tutorial]

- Application to Multiview Geometry:

  - part 1: [tensor diagrams][ap1]
  - part 2: (in construction)
  

Installation
------------

        $ sudo apt-get install haskell-platform libgsl0-dev liblapack-dev
        $ cabal update
        $ cabal install hTensor

Test
----

        $ ghci
        > import Numeric.LinearAlgebra.Exterior
        > printA "%4.0f" $ leviCivita 4 !"pqrs" * cov (leviCivita 4)!"qrsu"

        p^4 x u_4
                    u
            -6    0    0    0
        p    0   -6    0    0
             0    0   -6    0
             0    0    0   -6




[source]: http://hackage.haskell.org/package/hTensor
[tutorial]: http://dis.um.es/profesores/alberto/material/hTensor.pdf
[ap1]: http://dis.um.es/profesores/alberto/material/htmvg1.pdf

