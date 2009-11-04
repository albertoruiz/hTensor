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
    mlSolve, mlSolveHomog
) where

import Numeric.LinearAlgebra.Array
import Numeric.LinearAlgebra.Array.Util
import Numeric.LinearAlgebra.Array.Internal(mkNArray)
import Numeric.LinearAlgebra
import Data.List

-- import Debug.Trace
-- 
-- debug x = trace (show x) x
-- debug' m f x = trace (m ++ show (f x)) x


-- | Solution of the nonhomogenous linear system a x = b, where a and b are
-- general multidimensional arrays. The structure and dimension names
-- of the result are inferred from the arguments.
mlSolve :: Array Double -- ^ coefficients (a)
        -> Array Double -- ^ target       (b)
        -> Array Double -- ^ result       (x)
mlSolve a b = x where
    nx = names a \\ names b
    na = names a \\ nx
    nb = names b \\ names a
    aM = matrixator a na nx
    bM = matrixator b na nb
    xM = linearSolveSVD aM bM
    dx = filter ((`elem`(nx++nb)) .iName) (dims a ++ dims b)
    x = mkNArray dx (flatten xM)
    -- info = (rows aM, cols aM, cols bM)
    -- xM = debug' "system: " (const info) $ linearSolveSVD aM bM


-- | Solution of the homogenous linear system a x = 0, where a is a
-- general multidimensional array.
-- 
-- If the system is overconstrained we return the MSE solution. Otherwise we return a list of
-- (orthonormal) possible solutions.
mlSolveHomog :: Array Double   -- ^ coefficients (a)
             -> [Name]         -- ^ desired dimensions for the result
                               --   (a subset selected from the target, since
                               --   it makes no sense to have extra dimensions
                               --   in the resulting zero array).
             -> Int            -- ^ maximum number of solutions
             -> [Array Double] -- ^ basis for the solutions (x)
mlSolveHomog a nx' w = xs where
    nx = filter (`elem` nx') (names a)
    na = names a \\ nx
    aM = matrixator a na nx
    r = rows aM
    c = cols aM
    tc = ((r + w) `min` c) `max` (r+1)
    m = if r < c then fromRows $ take tc $ cycle $ toRows aM
                 else aM
    (_,_,v) = svd m   -- check rank?
    rd = if r >=c then c-1 -- overconstrained
                  else r
    xVs = take w $ drop rd (toColumns v)
    dx = filter ((`elem`(nx)).iName) (dims a)
    xs = map (mkNArray dx) xVs
--     info = (r, c, rows m, cols m)
--     xs = debug' "system: " (const info) $ map (mkNArray dx) xVs

