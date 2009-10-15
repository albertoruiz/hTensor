-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Packed.Array.Decomposition
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  provisional
--
-- Common multidimensional array decompositions. See the paper by Kolda & Balder.
--
-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Array.Decomposition (
    hosvd, truncateFactors
) where

import Numeric.LinearAlgebra.Array
import Numeric.LinearAlgebra.Array.Util
import Numeric.LinearAlgebra
import Data.List
import Numeric.LinearAlgebra.Array.Internal(firstIdx)
import Control.Parallel.Strategies

{- | Multilinear Singular Value Decomposition (or Tucker's method, see Lathauwer et al.).

    The first element in the result pair is a list with the core (head) and rotations so that
    t == product (fst (hsvd t)).

    The second element is a list of singular values along each mode, to give some idea about core structure.

    Rotations are (hopefully) computed in parallel.
-}
hosvd :: Array Double -> ([Array Double],[Vector Double])
hosvd t = (factors,ss)
    where factors = core!(map head dummies) : zipWith (!) (map (fromMatrix None None . trans) rs) axs
          (rs,ss) = unzip $ parMap rnf usOfSVD $ flats t
          n = length rs
          dummies = take n $ map return ['a'..'z'] \\ names t
          axs = zipWith (++) dummies (names t)
          core = product $ t!(map head dummies) : zipWith (!) (map (fromMatrix None None) rs) axs


-- get the matrices of the flattened tensor for all dimensions
flats t = map (snd . flip firstIdx t) (names t)


--check trans/ctrans
usOfSVD m = if rows m < cols m
        then let (s2,u) = eigSH' $ m <> ctrans m
                 s = sqrt (abs s2)
              in (u,s)
        else let (s2,v) = eigSH' $ ctrans m <> m
                 s = sqrt (abs s2)
                 u = m <> v <> pinv (diag s)
              in (u,s)


ttake ns t = (foldl1' (.) $ zipWith (onIndex.take) ns (names t)) t

-- | Truncate a decomposition from the desired number of principal components in each dimension.
truncateFactors :: [Int] -> [Array Double] -> [Array Double]
truncateFactors _ [] = []
truncateFactors ns (c:rs) = ttake ns c : zipWith f rs ns
    where f r n = onIndex (take n) (head (names r)) r
