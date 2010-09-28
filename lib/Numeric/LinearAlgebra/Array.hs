{-# LANGUAGE UndecidableInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.LinearAlgebra.Array
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  provisional
-- Portability :  portable
--
-- Simple multidimensional array with useful numeric instances.
--
-- Contractions only require equal dimension.
--
-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Array (
    None(..),
    Array,
    listArray,
    scalar,
    index,
    (!),(!>),(~>),
    (.*),
    printA
) where

import Numeric.LinearAlgebra.Array.Simple
import Numeric.LinearAlgebra.Array.Util
import Numeric.LinearAlgebra.Array.Internal(namesR)
import Numeric.LinearAlgebra.Array.Display(printA)
import Data.Packed(Vector)

-- | Create an 'Array' from a list of parts (@index = 'newIndex' 'None'@).
index :: Coord t => Name -> [Array t] -> Array t
index = newIndex None


-- | Element by element product.
infixl 7 .*
(.*) :: (Coord a, Compat i) => NArray i a -> NArray i a -> NArray i a
(.*) = zipArray (*)

instance (Coord t, Compat i) => Eq (NArray i t) where
    t1 == t2 = sameStructure t1 t2 && coords t1 == coords (reorder (namesR t1) t2)

instance (Show (NArray i t), Coord t, Compat i) => Num (NArray i t) where
    (+) = zipArray (+)
    (*) = (|*|)
    negate t = scalar (-1) * t
    fromInteger n = scalar (fromInteger n)
    abs _ = error "abs for arrays not defined"
    signum _ = error "signum for arrays not defined"

instance (Coord t, Compat i, Num (NArray i t)) => Fractional (NArray i t) where
    fromRational = scalar . fromRational
    (/) = zipArray (/)
    recip = mapArray recip

instance (Coord t, Compat i, Fractional (NArray i t), Floating t, Floating (Vector t)) => Floating (NArray i t) where
    sin   = mapArray sin
    cos   = mapArray cos
    tan   = mapArray tan
    asin  = mapArray asin
    acos  = mapArray acos
    atan  = mapArray atan
    sinh  = mapArray sinh
    cosh  = mapArray cosh
    tanh  = mapArray tanh
    asinh = mapArray asinh
    acosh = mapArray acosh
    atanh = mapArray atanh
    exp   = mapArray exp
    log   = mapArray log
    (**)  = zipArray (**)
    sqrt  = mapArray sqrt
    pi    = scalar pi
