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

    rename, (!),

    parts,
    newIndex,

    mapArray, zipArray, (|*|), smartProduct,

    extract, onIndex,

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
import Data.Packed(Container(..))
import Numeric.LinearAlgebra.Array.Simple
import Data.List(intersperse)

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

