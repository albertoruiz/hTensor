{-# LANGUAGE FlexibleContexts #-}
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
    -- * HOSVD
    hosvd, truncateFactors,
    -- * CP
    cpAuto, cpRun, cpInitRandom, cpInitSvd,
    -- * Utilities
    ALSParam(..), defaultParameters
) where

import Numeric.LinearAlgebra.Array
import Numeric.LinearAlgebra.Array.Util
import Numeric.LinearAlgebra.Array.Solve
import Numeric.LinearAlgebra hiding ((.*))
import Data.List
import System.Random
import Control.Parallel.Strategies

{- | Multilinear Singular Value Decomposition (or Tucker's method, see Lathauwer et al.).

    The first element in the result pair is a list with the core (head) and rotations so that
    t == product (fst (hsvd t)).

    The second element is a list of singular values along each mode,
    to give some idea about core structure.
-}
hosvd :: Array Double -> ([Array Double],[Vector Double])
hosvd t = (factors,ss)
    where factors = core!(map head dummies) : zipWith (!) (map (fromMatrix None None . trans) rs) axs
          (rs,ss) = unzip $ parMap rdeepseq usOfSVD $ flats t
          n = length rs
          dummies = take n $ map return ['a'..'z'] \\ names t
          axs = zipWith (++) dummies (names t)
          core = product $ t!(map head dummies) : zipWith (!) (map (fromMatrix None None) rs) axs


-- get the matrices of the flattened tensor for all dimensions
flats t = map (flip fibers t) (names t)


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

-- | Truncate a hosvd decomposition from the desired number of principal components in each dimension.
truncateFactors :: [Int] -> [Array Double] -> [Array Double]
truncateFactors _ [] = []
truncateFactors ns (c:rs) = ttake ns c : zipWith f rs ns
    where f r n = onIndex (take n) (head (names r)) r

------------------------------------------------------------------------

frobT = pnorm PNorm2 . coords

------------------------------------------------------------------------

unitRows [] = error "unitRows []"
unitRows (c:as) = foldl1' (.*) (c:xs) : as' where
    (xs,as') = unzip (map g as)
    g a = (x,a')
        where n = head (names a) -- hmmm
              rs = parts a n
              scs = map frobT rs
              x = diagT scs (order c) `rename` (names c)
              a' = (zipWith (.*) (map (scalar.recip) scs)) `onIndex` n $ a


{- | Basic CP optimization for a given rank. The result includes the obtained sequence of errors.

For example, a rank 3 approximation can be obtained as follows, where initialization
is based on the hosvd:

@
(y,errs) = cpRun (cpInitSvd (fst $ hosvd t) 3) 0.01 1E-6 t
@

-}
cpRun :: [Array Double] -- ^ starting point
      -> ALSParam None Double     -- ^ optimization parameters
      -> Array Double -- ^ input array
      -> ([Array Double], [Double]) -- ^ factors and error history
cpRun s0 params t = (unitRows $ head s0 : sol, errs) where
    (sol,errs) = mlSolve params [head s0] (tail s0) t



{- | Experimental implementation of the CP decomposition, based on alternating
     least squares. We try approximations of increasing rank, until the relative reconstruction error is below a desired percent of Frobenius norm (epsilon).

     The approximation of rank k is abandoned if the error does not decrease at least delta% in an iteration.

    Practical usage can be based on something like this:

@
cp finit delta epsilon t = cpAuto (finit t) delta epsilon t

cpS = cp (InitSvd . fst . hosvd)
cpR s = cp (cpInitRandom s)
@

     So we can write

@
 \-\- initialization based on hosvd
y = cpS 0.01 1E-6 t

 \-\- (pseudo)random initialization
z = cpR seed 0.1 0.1 t
@

-}
cpAuto :: (Int -> [Array Double]) -- ^ Initialization function for each rank
       -> ALSParam None Double    -- ^ optimization parameters
       -> Array Double -- ^ input array
       -> [Array Double] -- ^ factors
cpAuto finit params t = fst . head . filter ((<epsilon params). head . snd)
                      . map (\r->cpRun (finit r) params t) $ [1 ..]

----------------------

-- | cp initialization based on the hosvd
cpInitSvd :: [NArray None Double] -- ^ hosvd decomposition of the target array
          -> Int                  -- ^ rank
          -> [NArray None Double] -- ^ starting point
cpInitSvd (hos) k = d:as
    where c:rs = hos
          as = trunc (replicate (order c) k) rs
          d = diagT (replicate k 1) (order c) `rename` (names c)
          trunc ns xs = zipWith f xs ns
              where f r n = onIndex (take n . cycle) (head (names r)) r

cpInitSeq rs t k = ones:as where
    auxIndx = take (order t) $ map return ['a'..] \\ names t
    ones = diagT (replicate k 1) (order t) `rename` auxIndx
    ts = takes (map (*k) (sizes t)) rs
    as = zipWith4 f ts auxIndx (names t) (sizes t)
    f c n1 n2 p = (listArray [k,p] c) `rename` [n1,n2]

takes [] _ = []
takes (n:ns) xs = take n xs : takes ns (drop n xs)

-- | pseudorandom cp initialization from a given seed
cpInitRandom :: Int        -- ^ seed
             -> NArray i t -- ^ target array to decompose
             -> Int        -- ^ rank
             -> [NArray None Double] -- ^ random starting point
cpInitRandom seed = cpInitSeq (randomRs (-1,1) (mkStdGen seed))

----------------------------------------------------------------------
