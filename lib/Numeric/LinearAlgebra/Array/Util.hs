-- {-# LANGUAGE FlexibleInstances, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Packed.Array.Util
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  BSD3
-- Maintainer  :  Alberto Ruiz
-- Stability   :  provisional
--
-- Additional tools for manipulation of multidimensional arrays.
--
-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Array.Util (
    Coord, Compat(..),
    NArray, Idx(..), Name,
    scalar,
    order, names, size, sizes, typeOf, dims, coords,

    renameExplicit, (!>), renameO, (!),

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
    mkFun, mkAssoc, setType,
    renameParts,
    resetCoords,
    asScalar, asVector, asMatrix, applyAsMatrix,
    fibers, matrixator, matrixatorFree, analyzeProduct,
    fromVector, fromMatrix
    -- ,Container(..),
) where

import Numeric.LinearAlgebra.Array.Internal
import Numeric.LinearAlgebra.Array.Display
import Data.Packed(Matrix)
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
(!>) :: (Compat i, Coord t) => NArray i t -> [Char] -> NArray i t
infixl 9 !>
t !> s = renameExplicit (map f (words s)) t
  where
    f [a,b] = ([a],[b])
    f _ = error "impossible pattern in hTensor (!>)"

-- | Rename indices in alphabetical order. Equal indices of compatible type are contracted out.
renameO :: (Coord t, Compat i)
    => NArray i t
    -> [Name]
    -> NArray i t
renameO t ns = renameExplicit (zip od ns) t
    where od = map iName (sort (dims t))


-- | Rename indices in alphabetical order ('renameO') using single letter names.
(!) :: (Compat i, Coord t) => NArray i t -> [Char] -> NArray i t
infixl 9 !
t ! s = renameExplicit (zip od (map return s)) t
    where od = map iName (sort (dims t))


-- -- | 'renameRaw' the indices (in the internal order) with single-letter names. Equal indices of compatible type are contracted out.
-- infixl 8 !!!
-- (!!!) :: (Coord t, Compat i)
--        => NArray i t
--        -> String   -- ^ new indices
--        -> NArray i t
-- t !!! ns = renameRaw t (map return ns)


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

-- | Define an array using a function.
mkFun :: [Int] -> ([Int] -> Double) -> Array Double
mkFun ds f = listArray ds $ map f (sequence $ map (enumFromTo 0 . subtract 1. fromIntegral) $ ds)

-- | Define an array using an association list.
mkAssoc :: [Int] -> [([Int], Double)] -> Array Double
mkAssoc ds ps = mkFun ds f where
   f = maybe 0 id . flip lookup ps

-- | Change type of index.
setType :: (Compat i, Coord t) => Name -> i -> NArray i t -> NArray i t
setType n t a = mapDims f a where
    f i | iName i == n = i {iType = t}
        | otherwise    = i

-- | Extract the 'parts' of an array, and renameRaw one of the remaining indices
-- with succesive integers.
renameParts :: (Compat i, Coord t)
            => Name         -- ^ index of the parts to extract
            -> NArray i t   -- ^ input array
            -> Name         -- ^ index to renameRaw
            -> String       -- ^ prefix for the new names
            -> [NArray i t] -- ^ list or results
renameParts p t x pre = zipWith renameExplicit [[(x,pre ++ show k)] | k<-[1::Int ..] ] (parts t p)


applyAsMatrix :: (Coord t, Compat i) => (Matrix t -> Matrix t) -> (NArray i t -> NArray i t)
applyAsMatrix f t = flip renameRaw nms . fromMatrix r c . f . asMatrix $ t
    where [r,c] = map (flip typeOf t) nms
          nms = sort . namesR $ t
