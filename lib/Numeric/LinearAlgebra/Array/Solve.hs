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
    solve,
    solveHomog, solveHomog1, solveH,
    solveP,
-- *  Multilinear systems
-- ** General
    ALSParam(..), defaultParameters,
    mlSolve, mlSolveH,
-- ** Factorized
    solveFactors, solveFactorsH,
-- * Utilities
    eps, eqnorm, infoRank,
    solve', solveHomog', solveHomog1', solveP'
) where

import Numeric.LinearAlgebra.Array.Util
import Numeric.LinearAlgebra.Exterior
import Numeric.LinearAlgebra.Array.Internal(mkNArray, selDims, debug)
import Numeric.LinearAlgebra hiding ((.*))
import Data.List
import System.Random


-- | Solution of the linear system a x = b, where a and b are
-- general multidimensional arrays. The structure and dimension names
-- of the result are inferred from the arguments.
solve :: (Compat i, Coord t)
        => NArray i t -- ^ coefficients (a)
        -> NArray i t -- ^ target       (b)
        -> NArray i t -- ^ result       (x)
solve = solve' id

solve' g a b = x where
    nx = names a \\ names b
    na = names a \\ nx
    nb = names b \\ names a
    aM = g $ matrixator a na nx
    bM = matrixator b na nb
    xM = linearSolveSVD aM bM
    dx = map opos (selDims (dims a) nx) ++ selDims (dims b) nb
    x = mkNArray dx (flatten xM)


-- | Solution of the homogeneous linear system a x = 0, where a is a
-- general multidimensional array.
--
-- If the system is overconstrained we may provide the theoretical rank to get a MSE solution.
solveHomog :: (Compat i, Coord t)
           =>  NArray i t    -- ^ coefficients (a)
           -> [Name]         -- ^ desired dimensions for the result
                             --   (a subset selected from the target).
           -> Either Double Int -- ^ Left \"numeric zero\" (e.g. eps), Right \"theoretical\" rank
           -> [NArray i t] -- ^ basis for the solutions (x)
solveHomog = solveHomog' id

solveHomog' g a nx' hint = xs where
    nx = filter (`elem` (names a)) nx'
    na = names a \\ nx
    aM = g $ matrixator a na nx
    vs = nullspaceSVD hint aM (svd aM)
    dx = map opos (selDims (dims a) nx)
    xs = map (mkNArray dx) vs

-- | A simpler way to use 'solveHomog', which returns just one solution.
-- If the system is overconstrained it returns the MSE solution.
solveHomog1 :: (Compat i, Coord t)
            => NArray i t
            -> [Name]
            -> NArray i t
solveHomog1 = solveHomog1' id

solveHomog1' g m ns = head $ solveHomog' g m ns (Right (k-1))
    where k = product $ map iDim $ selDims (dims m) ns

-- | 'solveHomog1' for single letter index names.
solveH :: (Compat i, Coord t) => NArray i t -> [Char] -> NArray i t
solveH m ns = solveHomog1 m (map return ns)


-- | Solution of the linear system a x = b, where a and b are
-- general multidimensional arrays, with homogeneous equality along a given index.
solveP :: Tensor Double   -- ^ coefficients (a)
       -> Tensor Double   -- ^ desired result (b)
       -> Name            -- ^ the homogeneous dimension
       -> Tensor Double   -- ^ result (x)
solveP = solveP' id

solveP' g a b h = mapTat (solveP1 g h a) (names b \\ (h:names a)) b

-- solveP for a single right hand side
solveP1 g nh a b = solveHomog1' g ou ns where
    k = size nh b
    epsi = cov $ leviCivita k `rename` (nh : (take (k-1) $ (map (('e':).(:[])) ['2'..])))
    ou = a .* b' * epsi
    ns = (names a \\ names b) ++ x
    b' = renameExplicit [(nh,"e2")] b
    x = if nh `elem` (names a) then [] else [nh]

-----------------------------------------------------------------------

-- | optimization parameters for alternating least squares
data ALSParam i t = ALSParam
    { nMax  ::   Int     -- ^ maximum number of iterations
    , delta ::   Double  -- ^ minimum relative improvement in the optimization (percent, e.g. 0.1)
    , epsilon :: Double  -- ^ maximum relative error. For nonhomogeneous problems it is
                         --   the reconstruction error in percent (e.g.
                         --   1E-3), and for homogeneous problems is the frobenius norm of the
                         --  expected zero structure in the right hand side.
    , post :: [NArray i t] -> [NArray i t] -- ^ post-processing function after each full iteration (e.g. 'id')
    , postk :: Int -> NArray i t -> NArray i t-- ^ post-processing function for the k-th argument (e.g. 'const' 'id')
    , presys :: Matrix t -> Matrix t -- ^ preprocessing function for the linear systems (eg. 'id', or 'infoRank')
    }


optimize :: (x -> x)      -- ^ method
         -> (x -> Double) -- ^ error function
         -> x             -- ^ starting point
         -> ALSParam i t     -- ^ optimization parameters
         -> (x, [Double]) -- ^ solution and error history
optimize method errfun s0 p = (sol,e) where
    sols = take (max 1 (nMax p)) $ iterate method s0
    errs = map errfun sols
    (sol,e) = convergence (zip sols errs) []
    convergence [] _  = error "impossible"
    convergence [(s,err)] prev = (s, err:prev)
    convergence ((s1,e1):(s2,e2):ses) prev
        | e1 < epsilon p = (s1, e1:prev)
        | abs (100*(e1 - e2)/e1) < delta p = (s2, e2:prev)
        | otherwise = convergence ((s2,e2):ses) (e1:prev)

percent t s = 100 * frobT (t - smartProduct s) / frobT t

frobT t = pnorm PNorm2 . coords $ t

dropElemPos k xs = take k xs ++ drop (k+1) xs
replaceElemPos k v xs = take k xs ++ v : drop (k+1) xs

takes [] _ = []
takes (n:ns) xs = take n xs : takes ns (drop n xs)

-----------------------------------------------------------------------

-- | Solution of a multilinear system a x y z ... = b based on alternating least squares.
mlSolve
  :: (Compat i, Coord t, Num (NArray i t), Normed (Vector t))
     => ALSParam i t     -- ^ optimization parameters
     -> [NArray i t]  -- ^ coefficients (a), given as a list of factors.
     -> [NArray i t]  -- ^ initial solution [x,y,z...]
     -> NArray i t    -- ^ target (b)
     -> ([NArray i t], [Double]) -- ^ Solution and error history
mlSolve = als

als params a x0 b
    = optimize (post params .alsStep params a b) (percent b . (a++)) x0 params

alsStep params a b x = (foldl1' (.) (map (alsArg params a b) [0.. length x-1])) x

alsArg _ _ _ _ [] = error "alsArg _ _ []"
alsArg params a b k xs = sol where
    p = smartProduct (a ++ dropElemPos k xs)
    x = solve' (presys params) p b
    x' = postk params k x
    sol = replaceElemPos k x' xs

----------------------------------------------------------

-- | Solution of the homogeneous multilinear system a x y z ... = 0 based on alternating least squares.
mlSolveH
  :: (Compat i, Coord t, Num (NArray i t), Normed (Vector t))
     => ALSParam  i t    -- ^ optimization parameters
     -> [NArray i t]  -- ^ coefficients (a), given as a list of factors.
     -> [NArray i t]  -- ^ initial solution [x,y,z...]
     -> ([NArray i t], [Double]) -- ^ Solution and error history
mlSolveH = alsH

alsH params a x0
    = optimize (post params .alsStepH params a) (frobT . smartProduct . (a++)) x0 params

alsStepH params a x = (foldl1' (.) (map (alsArgH params a) [0.. length x-1])) x

alsArgH _ _ _ [] = error "alsArg _ _ []"
alsArgH params a k xs = sol where
    p = smartProduct (a ++ dropElemPos k xs)
    x = solveHomog1' (presys params) p (names (xs!!k))
    x' = postk params k x
    sol = replaceElemPos k x' xs

-------------------------------------------------------------

{- | Given two arrays a (source) and  b (target), we try to compute linear transformations x,y,z,... for each dimension, such that product [a,x,y,z,...] == b.
(We can use 'eqnorm' for 'post' processing, or 'id'.)
-}
solveFactors :: (Coord t, Random t, Compat i, Num (NArray i t), Normed (Vector t))
             => Int          -- ^ seed for random initialization
             -> ALSParam i t     -- ^ optimization parameters
             -> [NArray i t] -- ^ source (also factorized)
             -> String       -- ^ index pairs for the factors separated by spaces
             -> NArray i t   -- ^ target
             -> ([NArray i t],[Double]) -- ^ solution and error history
solveFactors seed params a pairs b =
    mlSolve params a (initFactorsRandom seed (smartProduct a) pairs b) b

initFactorsSeq rs a pairs b | ok = as
                            | otherwise = error "solveFactors index pairs"
  where
    (ia,ib) = unzip (map (\[x,y]->([x],[y])) (words pairs))
    ic = intersect (names a) (names b)
    ok = sort (names b\\ic) == sort ib && sort (names a\\ic) == sort ia
    db = selDims (dims b) ib
    da = selDims (dims a) ia
    nb = map iDim db
    na = map iDim da
    ts = takes (zipWith (*) nb na) rs
    as = zipWith5 f ts ib ia db da
    f c i1 i2 d1 d2 = (mkNArray [d1,opos d2] (fromList c)) `rename` [i1,i2]

initFactorsRandom seed a b = initFactorsSeq (randomRs (-1,1) (mkStdGen seed)) a b


-- | Homogeneous factorized system. Given an array a,
-- given as a list of factors as, and a list of pairs of indices
-- [\"pi\",\"qj\", \"rk\", etc.], we try to compute linear transformations
-- x!\"pi\", y!\"pi\", z!\"rk\", etc. such that product [a,x,y,z,...] == 0.
solveFactorsH
  :: (Coord t, Random t, Compat i, Num (NArray i t), Normed (Vector t))
     => Int -- ^ seed for random initialization
     -> ALSParam  i t    -- ^ optimization parameters
     -> [NArray i t] -- ^ coefficient array (a), (also factorized)
     -> String       -- ^ index pairs for the factors separated by spaces
     -> ([NArray i t], [Double]) -- ^ solution and error history
solveFactorsH seed params a pairs =
    mlSolveH params a (initFactorsHRandom seed (smartProduct a) pairs)

initFactorsHSeq rs a pairs = as where
    (ir,it) = unzip (map (\[x,y]->([x],[y])) (words pairs))
    nr = map (flip size a) ir
    nt = map (flip size a) it
    ts = takes (zipWith (*) nr nt) rs
    as = zipWith5 f ts ir it (selDims (dims a) ir) (selDims (dims a) it)
    f c i1 i2 d1 d2 = (mkNArray (map opos [d1,d2]) (fromList c)) `rename` [i1,i2]

initFactorsHRandom seed a pairs = initFactorsHSeq (randomRs (-1,1) (mkStdGen seed)) a pairs

----------------------------------

-- | post processing function that modifies a list of tensors so that they
-- have equal frobenius norm
eqnorm :: (Coord t, Coord (Complex t), Compat i, Num (NArray i t), Normed (Vector t) )
       => [NArray i t] -> [NArray i t]

eqnorm [] = error "eqnorm []"
eqnorm as = as' where
    n = length as
    fs = map (frobT) as
    s = product fs ** (1/fromIntegral n)
    as' = zipWith g as fs where g a f = a * real (scalar (s/f))

-- | nMax = 20, epsilon = 1E-3, delta = 1, post = id, postk = const id, presys = id
defaultParameters :: ALSParam i t
defaultParameters = ALSParam {
    nMax = 20,
    epsilon = 1E-3,
    delta = 1,
    post = id,
    postk = const id,
    presys = id
  }

-- | debugging function (e.g. for 'presys'), which shows rows, columns and rank of the
-- coefficient matrix of a linear system.
infoRank :: Field t => Matrix t -> Matrix t
infoRank a = debug "" (const (rows a, cols a, rank a)) a
