{-# LANGUAGE FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.LinearAlgebra.Exterior
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  experimental
--
-- Exterior Algebra.
--
--
-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Exterior (
    (/\),
    inner,
    leviCivita,
    dual,
    (\/),
    module Numeric.LinearAlgebra.Tensor,
    asMultivector, fromMultivector
) where

import Numeric.LinearAlgebra.Tensor
import Numeric.LinearAlgebra.Array.Internal
import Numeric.LinearAlgebra.Multivector(Multivector,fromTensor,maxDim,grade)
import qualified Numeric.LinearAlgebra.Multivector as MV
import Data.List

-- import Debug.Trace
-- debug x = trace (show x) x

interchanges :: (Ord a) => [a] -> Int
interchanges ls = sum (map (count ls) ls)
    where count l p = length $ filter (>p) $ take pel l
              where Just pel = elemIndex p l

signature :: (Num t, Ord a) => [a] -> t
signature l | length (nub l) < length l =  0
            | even (interchanges l)     =  1
            | otherwise                 = -1

gsym f t = mkNArray (dims t) (coords $ sum ts) where
    ns = map show [1 .. order t]
    t' = cov $ renameRaw t ns
    per = permutations ns
    ts  = map (flip renameRaw ns . f . flip reorder t') per

-- symmetrize t = gsym id t

antisymmetrize t = gsym scsig t
    where scsig x = scalar (signature (namesR x)) * x

fact n = product [1..n]

wedge a b = antisymmetrize (a*b) * (recip . fromIntegral) (fact (order a) * fact (order b))

infixl 5 /\
-- | The exterior (wedge) product of two tensors. Obtains the union of subspaces.
--
--   Implemented as the antisymmetrization of the tensor product.
(/\) :: (Coord t)
     => Tensor t
     -> Tensor t
     -> Tensor t
a /\ b = renseq (wedge a' b')
    where a' = renseq  a
          b' = renseq' b

-- levi n = antisymmetrize $ product $ zipWith renameRaw ts is
--     where is = map (return.show) [1 .. n]
--           ts = map (listTensor [n]) (toLists $ ident n)

levi n = listTensor (replicate n n) $ map signature $ sequence (replicate n [1..n])

-- | The full antisymmetric tensor of order n (contravariant version).
leviCivita :: Int -> Tensor Double
leviCivita = (map levi [0..] !!)

infixl 4 \/
-- | The \"meet\" operator. Obtains the intersection of subspaces.
--
-- @a \\\/ b = dual (dual a \/\\ dual b)@
(\/) :: Tensor Double -> Tensor Double -> Tensor Double
a \/ b = dual (dual a /\ dual b)

dual' n t = inner (leviCivita n) t

-- | Inner product of a r-vector with the whole space.
--
-- @dual t = inner (leviCivita n) t@
dual :: Tensor Double -> Tensor Double
dual t | isScalar t = error $ "cannot deduce dimension for dual of a scalar. Use s * leviCivita n"
       | otherwise  = dual' n t
    where n = case common iDim (dims t) of
                Just x -> x
                Nothing -> error $ "dual with different dimensions"

-- | Euclidean inner product of multivectors.
inner :: (Coord t)
      => Tensor t
      -> Tensor t
      -> Tensor t
inner a b | order a < order b = switch (renseq a) * renseq b * k
          | otherwise       = renseq a * switch (renseq b) * k
    where k = recip . fromIntegral $ fact $ min (order a) (order b)

renseq t = renameRaw t (map show [1..order t])
renseq' t = renameRaw t (map ((' ':).show) [1..order t])

isScalar = null . dims

-- | Extract a compact multivector representation from a full antisymmetric tensor.
--
-- asMultivector = Multivector.'fromTensor'.
--
-- (We do not check that the tensor is actually antisymmetric.)
asMultivector :: Tensor Double -> Multivector
asMultivector = fromTensor

-- | Create an explicit antisymmetric 'Tensor' from the components of a Multivector of a given grade.
fromMultivector :: Int -> Multivector -> Tensor Double
fromMultivector k t = sum $ map f (MV.coords $ grade k t) where
    f (x,es) = scalar x * foldl1' (/\) (map g es)
    n = maxDim t
    g i = vector $ replicate (i-1) 0 ++ 1 : replicate (n-i) 0
