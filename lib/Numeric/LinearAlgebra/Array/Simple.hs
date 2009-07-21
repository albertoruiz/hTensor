{-# LANGUAGE FlexibleInstances, FlexibleContexts, TypeSynonymInstances #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Packed.Array.Simple
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  provisional
-- Portability :  portable
--
-- Simple multidimensional arrays.
-- Contractions only require equal dimension.
--
-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Array.Simple (
    None(..),
    Array,
    listArray
) where

import Numeric.LinearAlgebra.Array.Internal
import Data.Packed


instance Show (Idx None) where
    show (Idx n s _t) = show n ++ ":" ++ s

-- | Unespecified coordinate type. Contractions only
-- require equal dimension.
data None = None deriving Eq


instance Compat None where
    compat d1 d2 = iDim d1 == iDim d2


-- | Multidimensional array with unespecified coordinate type.
type Array t = NArray None t

instance (Coord t) => Show (Array t) where
    show t | null (dims t) = "scalar "++ show (coords t @>0)
           | otherwise = "listArray "++ show (dims t) ++ " "++ show (toList $ coords t)

-- | Construction of an 'Array' from a list of dimensions and a list of elements in left to right order.
listArray :: (Coord t)
    => [Int] -- ^ dimensions
    -> [t]   -- ^ elements
    -> Array t
listArray ds cs = mkNArray dms (product ds |> (cs ++ repeat 0))
    where dms = zipWith3 Idx ds (map show [1::Int ..]) (repeat None)


