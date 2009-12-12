{-# LANGUAGE FlexibleInstances, FlexibleContexts, TypeSynonymInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.LinearAlgebra.Tensor
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  experimental
--
-- Tensor computations. Indices can only be contracted if they are of different 'Variant' type.
--
-----------------------------------------------------------------------------


module Numeric.LinearAlgebra.Tensor (
    -- * The Tensor type
    Tensor, Variant(..),
    listTensor,
    -- * Tensor creation utilities
    superindex, subindex,
    vector, covector, transf,
    -- * Index manipulation
    switch, cov, contrav, forget,
    -- * General array operations
    module Numeric.LinearAlgebra.Array
) where

import Numeric.LinearAlgebra.Array.Internal
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Array
import Data.List(intersperse)

type Tensor t = NArray Variant t

data Variant = Contra | Co deriving (Eq,Show)

instance Compat Variant where
    compat d1 d2 = iDim d1 == iDim d2 && iType d1 /= iType d2
    opos (Idx x n s) = Idx (flipV x) n s

instance Show (Idx Variant) where
    show (Idx Co n s)     = s ++ "_" ++ show n
    show (Idx Contra n s) = s ++ "^" ++ show n

instance (Coord t) => Show (Tensor t) where
    show t | null (dims t) = "scalar "++ show (coords t @>0)
           | order t == 1 = ixn ++ show n ++" " ++ (show . toList . coords $ t)
           | otherwise = ixn ++ show n ++ " [" ++ ps ++ "]"
      where n = head (names t)
            ps = concat $ intersperse ", " $ map show (parts t n)
            ixn = idxn (typeOf n t)
            idxn Co     = "subindex "
            idxn Contra = "superindex "


flipV Co = Contra
flipV Contra = Co

-- | Creates a tensor from a list of dimensions and a list of coordinates.
-- A positive dimension means that the index is assumed to be contravariant (vector-like), and
-- a negative dimension means that the index is assumed to be covariant (like a linear function, or covector). Contractions can only be performed between indices of different type.
listTensor :: Coord t
           => [Int] -- ^ dimensions
           -> [t]   -- ^ coordinates
           -> Tensor t
listTensor ds cs = mkNArray dms (product ds' |> (cs ++ repeat 0))
    where dms = zipWith3 Idx (map f ds) ds' (map show [1::Int ..])
          ds' = map abs ds
          f n | n>0       = Contra
              | otherwise = Co

-- | Create an 'Tensor' from a list of parts with a contravariant index (@superindex = 'newIndex' 'Contra'@).
superindex :: Coord t => Name -> [Tensor t] -> Tensor t
superindex = newIndex Contra

-- | Create an 'Tensor' from a list of parts with a covariant index (@subindex = 'newIndex' 'Co'@).
subindex :: Coord t => Name -> [Tensor t] -> Tensor t
subindex   = newIndex Co



-- | Change the 'Variant' nature of all dimensions to the opposite ones.
switch :: Tensor t -> Tensor t
switch = mapTypes flipV

-- | Make all dimensions covariant.
cov :: NArray i t -> Tensor t
cov     = mapTypes (const Co)

-- | Make all dimensions contravariant.
contrav :: NArray i t -> Tensor t
contrav = mapTypes (const Contra)

-- | Remove the 'Variant' nature of coordinates.
forget :: NArray i t -> Array t
forget = mapTypes (const None)

--------------------------------------------------------------

-- | Create a contravariant 1st order tensor from a list of coordinates.
vector :: [Double] -> Tensor Double
vector   = fromVector Contra . fromList

-- | Create a covariant 1st order tensor from a list of coordinates.
covector :: [Double] -> Tensor Double
covector = fromVector Co . fromList

-- | Create a 1-contravariant, 1-covariant 2nd order from list of lists of coordinates.
transf :: [[Double]] -> Tensor Double
transf   = fromMatrix Contra Co . fromLists
