-- {-# LANGUAGE FlexibleInstances, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Packed.Array.Util
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  provisional
-- Portability :  portable
--
-- Additional tools for manipulation of multidimensional arrays.
--
-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Array.Util (
    Coord, Compat(..),
    NArray, Idx(..), Name,
    scalar,
    order, names, size, sizes, typeOf, dims, coords,

    rename, (!), renameExplicit, (>@>), (!>),

    parts,
    newIndex,

    mapArray, zipArray, (|*|), smartProduct, outers,

    extract, onIndex, mapTat,

    reorder, (~>),
    formatArray, formatFixed, formatScaled,
    dummyAt, noIdx,
    conformable,
    sameStructure,
    makeConformant,
    basisOf,
    atT, takeDiagT, diagT,
    asScalar, asVector, asMatrix,
    fibers, matrixator, analyzeProduct,
    fromVector, fromMatrix,
    Container(..),
) where

import Numeric.LinearAlgebra.Array.Internal
import Numeric.LinearAlgebra.Array.Display
import Data.Packed(Container(..))
import Numeric.LinearAlgebra.Array.Simple
import Data.List(intersperse,sort,foldl1')

-- infixl 9 #
-- (#) :: [Int] -> [Double] -> Array Double
-- (#) = listArray

-- | Multidimensional diagonal of given order.
diagT :: [Double] -> Int -> Array Double
diagT v n = replicate n k `listArray` concat (intersperse z (map return v))
    where k = length v
          tot = k^n
          nzeros = (tot - k) `div` (k-1)
          z = replicate nzeros 0


-- | Explicit renaming of single letter index names.
--
-- For instance, @t >\@> \"pi qj\"@ changes index \"p\" to \"i\" and \"q\" to \"j\".
(>@>) :: (Compat i, Coord t) => NArray i t -> [Char] -> NArray i t
infixl 9 >@>
t >@> s = renameExplicit (map (\[a,b]->([a],[b])) (words s)) t


-- | Rename the ordered dimensions of an array (contravariant < covariant), using single letter names.
(!>) :: (Compat i, Coord t) => NArray i t -> [Char] -> NArray i t
infixl 9 !>
t !> s = renameExplicit (zip od (map return s)) t
    where od = map iName (sort (dims t))


-- | 'rename' the indices (in the internal order) with single-letter names. Equal indices of compatible type are contracted out.
infixl 8 !
(!) :: (Coord t, Compat i)
       => NArray i t
       -> String   -- ^ new indices
       -> NArray i t
t ! ns = rename t (map return ns)


-- | 'reorder' (transpose) dimensions of the array (with single letter names).
--
-- Operations are defined by named indices, so the transposed array is operationally equivalent to the original one.
infixl 8 ~>
(~>) :: (Coord t) => NArray i t -> String -> NArray i t
t ~> ns = reorder (map return ns) t


-- | Map a function at the internal level selected by a set of indices
mapTat :: (Coord a, Coord b, Compat i)
         => (NArray i a -> NArray i b)
         -> [Name]
         -> NArray i a
         -> NArray i b
mapTat f [] = f
mapTat f (a:as) = onIndex (map $ mapTat f as) a

-- | Outer product of a list of arrays along the common indices.
outers :: (Coord a, Compat i) => [NArray i a] -> NArray i a
outers = foldl1' (zipArray (*))
