-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Packed.Array.Solve
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  provisional
--
-- Solution of a general multidimensional linear system a x = b.
--
-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Array.Solve (
    mlSolve, mlSolveHomog, mlSolveH,
    als
) where

import Numeric.LinearAlgebra.Array.Util
import Numeric.LinearAlgebra.Array.Internal(mkNArray)
import Numeric.LinearAlgebra
import Data.List

-- import Debug.Trace
-- 
-- debug x = trace (show x) x
-- debug' m f x = trace (m ++ show (f x)) x
-- 
-- sv m = let (_,s,_) = svd m in s


-- | Solution of the nonhomogenous linear system a x = b, where a and b are
-- general multidimensional arrays. The structure and dimension names
-- of the result are inferred from the arguments.
mlSolve :: (Compat i, Coord t)
        => NArray i t -- ^ coefficients (a)
        -> NArray i t -- ^ target       (b)
        -> NArray i t -- ^ result       (x)
mlSolve a b = x where
    nx = names a \\ names b
    na = names a \\ nx
    nb = names b \\ names a
    aM = {-debug' "aMnew: " id $-} matrixator a na nx
    bM = {-debug' "bMnew: " id $-} matrixator b na nb
    xM = linearSolveSVD ({-debug' "info = " (const info)-} aM) bM
    dx = map opos (selDims (dims a) nx) ++ selDims (dims b) nb
    -- filter ((`elem`(nx++nb)) .iName) (dims a ++ dims b)
    x = mkNArray dx (flatten xM)
--    info = ((rows aM, cols aM, cols bM), (na, nx, nb), dims a, dims b)
--    xM = debug' "system: " (const info) $ linearSolveSVDR (Just (1E3*eps)) aM bM
--    info = show (names a) ++ show (names b) ++ show na ++ show nx ++ show nb ++ show dx

selDims ds = map f where
    f n = head $ filter ((n==).iName) ds

-- | Solution of the homogeneous linear system a x = 0, where a is a
-- general multidimensional array.
--
-- If the system is overconstrained we may provide the theoretical rank to get a MSE solution.
mlSolveHomog :: (Compat i, Coord t)
             => NArray i t     -- ^ coefficients (a)
             -> [Name]         -- ^ desired dimensions for the result
                               --   (a subset selected from the target, since
                               --   it makes no sense to have extra dimensions
                               --   in the resulting zero array).
             -> Either Double Int -- ^ Left \"numeric zero\" (e.g. eps), Right \"theoretical\" rank
             -> [NArray i t] -- ^ basis for the solutions (x)
mlSolveHomog a nx' hint = xs where
    nx = filter (`elem` (names a)) nx'
    na = names a \\ nx
    aM = matrixator a na nx
    vs = nullspaceSVD hint aM (svd aM)
    dx = map opos (selDims (dims a) nx)
    xs = map (mkNArray dx) vs
    -- debug' "mlSolveHomog: " (const (rows aM, cols aM, rank aM)) $


mlSolveH m ns' = head $ mlSolveHomog m ns (Right (k-1))
    where k = product $ map iDim $ selDims (dims m) ns
          ns = map return ns'

-- mlSolveHomog' a nx' w = xs where
--     nx = filter (`elem` (names a)) nx'
--     na = names a \\ nx
--     aM = matrixator a na nx
--     r = rows aM
--     c = cols aM
--     k = rank aM
--     nsol = c - k `max` 1
--     tc = ((r + w) `min` c) `max` (r+1)
--     m = if r < c then fromRows $ take tc $ cycle $ toRows aM
--                  else aM
--     (_,_,v) = (debug' "mlSolveHomog: " (const (r,c,k,rcond aM))) svd m   -- check rank?
--     rd = if k ==c then c-1 -- overconstrained
--                   else k
--     xVs = take w $ drop rd (toColumns v)
--     dx = selDims (dims a) nx
--     xs = map (mkNArray dx) xVs
-- --     info = (r, c, rows m, cols m)
-- --     xs = debug' "system: " (const info) $ map (mkNArray dx) xVs


-----------------------------------------------------------------------

optimize :: (x -> x)      -- ^ method
         -> (x -> Double) -- ^ error function
         -> x             -- ^ starting point
         -> Double        -- ^ delta: minimum relative improvement in the optimization (percent, e.g. 0.1)
         -> Double -- ^ epsilon: desired error
         -> (x, [Double]) -- ^ solution and error history
optimize method errfun s0 delta epsilon = (sol,e) where
    sols = iterate method s0
    errs = map errfun sols
    (sol,e) = convergence (zip sols errs) []
    convergence [] _  = error "impossible"  -- to avoid warning
    convergence [_] _ = error "impossible"
    convergence ((s1,e1):(s2,e2):ses) prev
        | e1 < epsilon = (s1, e1:prev)
        | abs (100*(e1 - e2)/e1) < delta = (s2, e2:prev)
        | otherwise = convergence ((s2,e2):ses) (e1:prev)

percent t s = 100 * frobT (t - product s) / frobT t

frobT t = pnorm PNorm2 . coords $ t

-----------------------------------------------------------------------

als post t s0 delta epsilon
    = optimize (post.alsStep t) (percent t) s0 delta epsilon

alsStep t s = (foldl1' (.) (map (alsArg t) [1..n])) s
    where n = length s - 1

alsArg _ _ [] = error "alsArg _ _ []"
alsArg t k (d:as) = sol where
    p = product (d : dropElemPos k as)
    x = mlSolve p ({-debug' ": " dims $-} t)
    sol = d: replaceElemPos k x as

dropElemPos k xs = take (k-1) xs ++ drop k xs
replaceElemPos k v xs = take (k-1) xs ++ v : drop k xs
