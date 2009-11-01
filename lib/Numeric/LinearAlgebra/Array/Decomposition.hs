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
    diagT, takeDiagT
) where

import Numeric.LinearAlgebra.Array
import Numeric.LinearAlgebra.Array.Util
import Numeric.LinearAlgebra hiding ((.*))
import Numeric.LinearAlgebra.LAPACK (linearSolveLSR)
import Data.List
import Control.Parallel.Strategies
import System.Random

-- import Debug.Trace
-- 
-- debug x = trace (show x) x

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

cpArg _ _ [] = error "cpArg _ _ []"
cpArg t k (d:as) = d': replaceElemPos k ar as where
    [n1,n2] = names (as!!(k-1))
    namesta = replaceElem n2 n1 (names t)
    ta = reorder namesta (product $ d : dropElemPos k as)
    fa = trans $ fibers n1 ta
    ft = trans $ fibers n2 t
    a = toRows $ linearSolveLSR fa ft
    scs = map (pnorm PNorm2) a
    d' = d .* (diagT scs (order t) `rename` (names d))
    au = fromRows $ zipWith (*/) a scs
    ar = fromMatrix None None au `rename` [n1,n2]

-- | Basic optimization step for the CP approximation
cpStep :: Array Double -> [Array Double] -> [Array Double]
cpStep t = foldl1' (.) (map (cpArg t) [1..order t])

convergence _ _ [] _ = error "convergence on finite sequence"
convergence _ _ [_] _ = error "convergence on finite sequence"
convergence epsrel epsabs ((s1,e1):(s2,e2):ses) prev
    | e1 < epsabs = (s1, e1:prev)
    | abs (100*(e1 - e2)/e1) < epsrel = (s2, e2:prev)
    | otherwise = convergence epsrel epsabs ((s2,e2):ses) (e1:prev)

sizes = map iDim . dims

{- | Basic CP optimization for a given rank. The result includes the obtained sequence of errors.

For example, a rank 3 approximation can be obtained as follows, where initialization
is based on the hosvd:

@
(y,errs) = cpRun (cpInitSvd (fst $ hosvd t) 3) 0.01 1E-6 t
@

-}
cpRun :: [Array Double] -- ^ starting point
      -> Double -- ^ delta: minimum relative improvement in the optimization (percent, e.g. 0.1)
      -> Double -- ^ epsilon: desired relative reconstruction error (percent, e.g. 0.1)
      -> Array Double -- ^ input array
      -> ([Array Double], [Double]) -- ^ factors and error history
cpRun s0 delta epsilon t = (sol,e) where
    sols = iterate (cpStep t) s0
    errs = map (\s -> 100 * frobT (t - product s) / frobT t) sols
    (sol,e) = convergence delta epsilon (zip sols errs) []

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
       -> Double -- ^ delta: minimum relative improvement in the optimization (percent, e.g. 0.1)
       -> Double -- ^ epsilon: desired relative reconstruction error (percent, e.g. 0.1)
       -> Array Double -- ^ input array
       -> [Array Double] -- ^ factors
cpAuto finit delta epsilon t = fst . head . filter ((<epsilon). head . snd)
                             . map (\r->cpRun (finit r) delta epsilon t) $ [1 ..]

----------------------

dropElemPos k xs = take (k-1) xs ++ drop k xs

replaceElemPos k v xs = take (k-1) xs ++ v : drop k xs

replaceElem x v xs = map f xs where
    f a | x==a      = v
        | otherwise = a

frobT = pnorm PNorm2 . coords

infixl 9 #
(#) :: [Int] -> [Double] -> Array Double
(#) = listArray

-- | Multidimensional diagonal of given order.
diagT :: [Double] -> Int -> Array Double
diagT v n = replicate n k # concat (intersperse z (map return v))
    where k = length v
          tot = k^n
          nzeros = (tot - k) `div` (k-1)
          z = replicate nzeros 0

takeDiagT :: Coord t => Array t -> [t]
takeDiagT t = map (asScalar . atT t) cds where
    n = minimum (sizes t)
    o = order t
    cds = map (replicate o) [0..n-1]

atT :: Coord t => Array t -> [Int] -> Array t
atT t c = atT' c t where
    atT' cs = foldl1' (.) (map fpart cs)
    fpart k q = parts q (head (names q)) !! k

-- dropRank [] = []
-- dropRank (h:xs) = h':xs' where
--     k = posMin (map abs $ takeDiagT h)
--     h' = (foldl1' (.) $ map (\n-> dropElemPos k `onIndex` n) (names h)) h
--     xs' = map f xs
--         where f x = onIndex (dropElemPos k) (head (names x)) x
--     posMin ys = p+1 where Just p = elemIndex (minimum ys) ys

-- | cp inicialization based on the hosvd
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
    f c n1 n2 p = ([k,p] # c) `rename` [n1,n2]
    takes [] _ = []
    takes (n:ns) xs = take n xs : takes ns (drop n xs)

-- | pseudorandom cp inicialization from a given seed
cpInitRandom :: Int        -- ^ seed
             -> NArray i t -- ^ target array to decompose
             -> Int        -- ^ rank
             -> [NArray None Double] -- ^ random starting point
cpInitRandom seed = cpInitSeq (randomRs (-1,1) (mkStdGen seed))
