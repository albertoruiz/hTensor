-- {-# LANGUAGE FlexibleInstances, FlexibleContexts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Packed.Array.Util
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  provisional
-- Portability :  portable
--
-- Additional tools for manipulation of multidimensional arrays.
--
-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Array.Util (
    Coord, Compat(..),
    NArray, Idx(..), Name,
    scalar,
    rank, names, size, typeOf, dims, coords,

    rename, (!),

    parts,
    newIndex,

    mapArray, zipArray, (|*|),

    extract, onIndex,

    reorder, (~>),
    formatArray, formatFixed, formatScaled,
    dummyAt, noIdx,
    conformable,
    sameStructure,
    makeConformant,
    basisOf,
    asScalar, asVector, asMatrix,
    fromVector, fromMatrix,
    Container(..),
) where

import Numeric.LinearAlgebra.Array.Internal
import Data.Packed(Container(..))
