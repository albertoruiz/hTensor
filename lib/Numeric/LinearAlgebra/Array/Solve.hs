{-# LANGUAGE FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Packed.Array.Solve
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  provisional
--
-- Solution of general multidimensional linear and multilinear systems.
--
-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Array.Solve (
-- * Linear systems
    solve, solveHomog, solveH,
-- *  Multilinear systems
    mlSolve, mlSolveH,
-- * Utilities
    eps
) where

import Numeric.LinearAlgebra.Array.Util
import Numeric.LinearAlgebra.Array.Internal(mkNArray, selDims, debug)
import Numeric.LinearAlgebra
import Data.List


-- | Solution of the linear system a x = b, where a and b are
-- general multidimensional arrays. The structure and dimension names
-- of the result are inferred from the arguments.
solve :: (Compat i, Coord t)
        => NArray i t -- ^ coefficients (a)
        -> NArray i t -- ^ target       (b)
        -> NArray i t -- ^ result       (x)
solve a b = x where
    nx = names a \\ names b
    na = names a \\ nx
    nb = names b \\ names a
    aM = {-debug "aMnew: " id $-} matrixator a na nx
    bM = {-debug "bMnew: " id $-} matrixator b na nb
    xM = linearSolveSVD ({-debug "info = " (const info)-} aM) bM
    dx = map opos (selDims (dims a) nx) ++ selDims (dims b) nb
    -- filter ((`elem`(nx++nb)) .iName) (dims a ++ dims b)
    x = mkNArray dx (flatten xM)
    -- info = ((rows aM, cols aM, cols bM), (na, nx, nb), dims a, dims b)
--    xM = debug "system: " (const info) $ linearSolveSVDR (Just (1E3*eps)) aM bM
--    info = show (names a) ++ show (names b) ++ show na ++ show nx ++ show nb ++ show dx



-- | Solution of the homogeneous linear system a x = 0, where a is a
-- general multidimensional array.
--
-- If the system is overconstrained we may provide the theoretical rank to get a MSE solution.
solveHomog :: (Compat i, Coord t)
           => NArray i t     -- ^ coefficients (a)
           -> [Name]         -- ^ desired dimensions for the result
                             --   (a subset selected from the target, since
                             --   it makes no sense to have extra dimensions
                             --   in the resulting zero array).
           -> Either Double Int -- ^ Left \"numeric zero\" (e.g. eps), Right \"theoretical\" rank
           -> [NArray i t] -- ^ basis for the solutions (x)
solveHomog a nx' hint = xs where
    nx = filter (`elem` (names a)) nx'
    na = names a \\ nx
    aM = matrixator a na nx
    vs = nullspaceSVD hint aM (svd aM)
    dx = map opos (selDims (dims a) nx)
    xs =
         debug "mlSolveHomog: " (const (rows aM, cols aM, rank aM)) $
         map (mkNArray dx) vs

-- | A simpler way to use 'solveHomog' for single letter index names, which returns just one solution.
-- If the system is overconstrained it returns the MSE solution.
solveH :: (Compat i, Coord t) => NArray i t -> [Char] -> NArray i t
solveH m ns = solveH' m (map return ns)

solveH' m ns = head $ solveHomog m ns (Right (k-1))
    where k = product $ map iDim $ selDims (dims m) ns


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
--     (_,_,v) = (debug "mlSolveHomog: " (const (r,c,k,rcond aM))) svd m   -- check rank?
--     rd = if k ==c then c-1 -- overconstrained
--                   else k
--     xVs = take w $ drop rd (toColumns v)
--     dx = selDims (dims a) nx
--     xs = map (mkNArray dx) xVs
-- --     info = (r, c, rows m, cols m)
-- --     xs = debug "system: " (const info) $ map (mkNArray dx) xVs


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

dropElemPos k xs = take k xs ++ drop (k+1) xs
replaceElemPos k v xs = take k xs ++ v : drop (k+1) xs

-----------------------------------------------------------------------

-- | Solution of a multilinear system a x y z ... = b based on alternating least squares.
mlSolve
  :: (Compat i, Coord t, Num (NArray i t), Normed (Vector t)) =>
     ([NArray i t] -> [NArray i t])  -- ^ post-processing function after each iteration (e.g. id)
     -> Double        -- ^ delta: minimum relative improvement in the optimization (percent, e.g. 0.1)
     -> Double        -- ^ epsilon: desired relative reconstruction error (percent, e.g. 1E-3)
     -> NArray i t    -- ^ coefficients (a)
     -> [NArray i t]  -- ^ initial solution [x,y,z...]
     -> NArray i t    -- ^ target (b)
     -> ([NArray i t], [Double]) -- ^ Solution and error history
mlSolve = als

als post delta epsilon a x0 b
    = optimize (post.alsStep a b) (percent b . (a:)) x0 delta epsilon

alsStep a b x = (foldl1' (.) (map (alsArg a b) [0.. length x-1])) x

alsArg _ _ _ [] = error "alsArg _ _ []"
alsArg a b k xs = sol where
    p = product (a : dropElemPos k xs)
    x = solve p b
    sol = replaceElemPos k x xs

----------------------------------------------------------

-- | Solution of the homogeneous multilinear system a x y z ... = 0 based on alternating least squares.
mlSolveH
  :: (Compat i, Coord t, Num (NArray i t), Normed (Vector t)) =>
     ([NArray i t] -> [NArray i t])  -- ^ post-processing function after each iteration (e.g. id)
     -> Double        -- ^ delta: minimum relative improvement in the optimization (percent, e.g. 0.1)
     -> Double        -- ^ epsilon: frob norm of the right hand side
     -> NArray i t    -- ^ coefficients (a)
     -> [NArray i t]  -- ^ initial solution [x,y,z...]
     -> ([NArray i t], [Double]) -- ^ Solution and error history
mlSolveH = alsH

alsH post delta epsilon a x0
    = optimize (post.alsStepH a) (frobT . product . (a:)) x0 delta epsilon

alsStepH a x = (foldl1' (.) (map (alsArgH a) [0.. length x-1])) x

alsArgH _ _ [] = error "alsArg _ _ []"
alsArgH a k xs = sol where
    p = product (a : dropElemPos k xs)
    x = solveH' p (names (xs!!k))
    sol = replaceElemPos k x xs
