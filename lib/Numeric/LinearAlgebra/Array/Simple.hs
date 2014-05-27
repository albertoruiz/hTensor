{-# LANGUAGE FlexibleInstances, FlexibleContexts, TypeSynonymInstances #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Packed.Array.Simple
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  BSD3
-- Maintainer  :  Alberto Ruiz
-- Stability   :  provisional
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
import Data.List(intersperse)


instance Show (Idx None) where
    show (Idx _t n s) = s ++ ":" ++ show n

-- | Unespecified coordinate type. Contractions only
-- require equal dimension.
data None = None deriving (Eq,Show)


instance Compat None where
    compat d1 d2 = iDim d1 == iDim d2
    opos = id


-- | Multidimensional array with unespecified coordinate type.
type Array t = NArray None t

instance (Coord t) => Show (Array t) where
    show t | null (dims t) = "scalar "++ show (coords t @>0)
           | order t == 1 = "index " ++ show n ++" " ++ (show . toList . coords $ t)
           | otherwise = "index "++ show n ++ " [" ++ ps ++ "]"
      where n = head (namesR t)
            ps = concat $ intersperse ", " $ map show (parts t n)

-- ++ " "++ show (toList $ coords t)

-- | Construction of an 'Array' from a list of dimensions and a list of elements in left to right order.
listArray :: (Coord t)
    => [Int] -- ^ dimensions
    -> [t]   -- ^ elements
    -> Array t
listArray ds cs = mkNArray dms (product ds |> (cs ++ repeat 0))
    where dms = zipWith3 Idx (repeat None) ds (map show [1::Int ..])


