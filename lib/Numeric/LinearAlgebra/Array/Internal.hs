{-# OPTIONS_HADDOCK hide #-}
{-# LANGUAGE FlexibleInstances, FlexibleContexts, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Packed.Array.Internal
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  provisional
-- Portability :  portable
--
-- Multidimensional arrays.
--
-- The arrays provided by this library are immutable, built on top of hmatrix
-- structures.
-- Operations work on complete structures (indexless), and dimensions have \"names\", 
-- in order to select the desired contractions in tensor computations.
--
-- This module contains auxiliary functions not required by the end user.

-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Array.Internal (
    -- * Data structures
    NArray, Idx(..), Name,
    order, namesR, names, size, sizesR, sizes, typeOf , dims, coords,
    Compat(..),
    -- * Array creation
    scalar,
    mkNArray,
    fromVector, fromMatrix,
    -- * Array manipulation
    renameRaw,
    parts, partsRaw,
    (|*|),
    analyzeProduct,
    smartProduct,
    zipArray,
    mapArray,
    extract,
    onIndex,
    -- * Utilities
    seqIdx,
    reorder,
    sameStructure,
    conformable,
    makeConformant,
    mapTypes, mapNames,
    renameSuperRaw, renameExplicit,
    newIndex,
    basisOf,
    common,
    selDims, mapDims,
    takeDiagT, atT,
    firstIdx, fibers, matrixator, matrixatorFree,
    Coord,
    asMatrix, asVector, asScalar,
    debug
) where

import Data.Packed
import Data.List
import Numeric.LinearAlgebra((<>),Field)
import Control.Applicative
import Data.Function(on)
-- import Control.Parallel.Strategies
import Debug.Trace

debug m f x = trace (m ++ show (f x)) x

-- | Types that can be elements of the multidimensional arrays.
class (Num (Vector t), Field t) => Coord t
instance Coord Double
instance Coord (Complex Double)

-- | indices are denoted by strings, (frequently single-letter)
type Name = String

-- | Dimension descriptor.
data Idx i = Idx { iType :: i
                 , iDim  :: Int
                 , iName :: Name
                 } deriving (Eq)

instance Eq i => Ord (Idx i) where
    compare = compare `on` iName

-- | A multidimensional array with index type i and elements t.
data NArray i t = A { dims   :: [Idx i]   -- ^ Get detailed dimension information about the array.
                    , coords :: Vector t  -- ^ Get the coordinates of an array as a
                                          -- flattened structure (in the order specified by 'dims').
                    }

-- | development function not intended for the end user
mkNArray :: [Idx i] -> Vector a -> NArray i a
mkNArray [] _ = error "array with empty dimensions, use scalar"
mkNArray dms vec = A dms v where
    ds = map iDim dms
    n = product ds
    v = if dim vec == n && minimum ds > 0
            then vec
            else error $ show ds ++ " dimensions and " ++
                         show (dim vec) ++ " coordinates for mkNArray"

-- | Create a 0-dimensional structure.
scalar :: Coord t => t -> NArray i t
scalar x = A [] (fromList [x])

-- | Rename indices (in the internal order). Equal indices are contracted out.
renameRaw :: (Coord t, Compat i)
       => NArray i t
       -> [Name]     -- ^ new names
       -> NArray i t
renameRaw t ns = contract (renameSuperRaw t ns)

renameSuperRaw (A d v) l
    | length l == length d = A d' v
    | otherwise = error $ "renameRaw " ++ show d ++ " with " ++ show l
   where d' = zipWith f d l
         f i n = i {iName=n}

mapDims f (A d v) = A (map f d) v

mapTypes :: (i1 -> i2) -> NArray i1 t -> NArray i2 t
mapTypes f = mapDims (\i -> i {iType = f (iType i)})

mapNames :: (Name -> Name) -> NArray i t -> NArray i t
mapNames f = mapDims (\i -> i {iName = f (iName i)})

-- | Rename indices using an association list.
renameExplicit :: (Compat i, Coord t) => [(Name,Name)] -> NArray i t -> NArray i t
renameExplicit al = g . mapNames f where
    f n = maybe n id (lookup n al)
    g t = reorder orig (contract t) where orig = nub (namesR t) \\ common1 t

-- | Index names (in internal order).
namesR :: NArray i t -> [Name]
namesR = map iName . dims

-- | Index names.
names :: NArray i t -> [Name]
names = sort . namesR

-- | Dimension of given index.
size :: Name -> NArray i t -> Int
size n t = (iDim . head) (filter ((n==).iName) (dims t))

sizesR :: NArray i t -> [Int]
sizesR = map iDim . dims

sizes :: NArray i t -> [Int]
sizes t = map (flip size t) (names t)

-- | Type of given index.
typeOf :: Compat i => Name -> NArray i t -> i
typeOf n t = (iType . head) (filter ((n==).iName) (dims t))

-- | The number of dimensions of a multidimensional array.
order :: NArray i t -> Int
order = length . dims

selDims ds = map f where
    f n = head $ filter ((n==).iName) ds

----------------------------------------------------------

common2 t1 t2 = [ n1 | n1 <- namesR t1, n2 <- namesR t2, n1==n2]

analyzeProduct :: (Coord t, Compat i) => NArray i t -> NArray i t -> Maybe (NArray i t, Int)
analyzeProduct a b = r where
    nx  = common2 a b
    dx1 = selDims (dims a) nx
    dx2 = selDims (dims b) nx
    ok  = and $ zipWith compat dx1 dx2
    (tma,na) = matrixatorFree a nx
    (mb,nb)  = matrixatorFree b nx
    mc  = trans tma <> mb
    da  = selDims (dims a) na
    db  = selDims (dims b) nb
    dc  = da ++ db
    c   = A dc (flatten mc)
    sz  = product (map iDim dc)
    r | ok = Just (c, sz)
      | otherwise = Nothing

infixl 5 |*|
-- | Tensor product with automatic contraction of repeated indices, following Einstein summation convention.
(|*|) :: (Coord t, Compat i) => NArray i t -> NArray i t -> NArray i t
t1 |*| t2 = case analyzeProduct t1 t2 of
    Nothing -> error $ "wrong contraction2: "++(show $ dims t1)++" and "++(show $ dims t2)
    Just (r,_) -> r

----------------------------------------------------------

lastIdx name t = ((d1,d2),m) where
    (d1,d2) = span (\d -> iName d /= name) (dims t)
    c = product (map iDim d2)
    m = reshape c (coords t)

firstIdx name t = (nd,m')
    where ((d1,d2),m) = lastIdx name t
          m' = reshape c $ flatten $ trans m
          nd = d2++d1
          c = dim (coords t) `div` (iDim $ head d2)

-- | Obtain a matrix whose columns are the fibers of the array in the given dimension. The column order depends on the selected index (see 'matrixator').
fibers :: Coord t => Name -> NArray i t -> Matrix t
fibers n = snd . firstIdx n

-- | Reshapes an array as a matrix with the desired dimensions as flattened rows and flattened columns.
matrixator :: (Coord t) => NArray i t -- ^ input array
                        -> [Name]    -- ^ row dimensions
                        -> [Name]    -- ^ column dimensions
                        -> Matrix t   -- ^ result
matrixator t nr nc = reshape s (coords q) where
    q = reorder (nr++nc) t
    s = product (map (flip size t) nc)

-- | Reshapes an array as a matrix with the desired dimensions as flattened rows and flattened columns. We do not force the order of the columns.
matrixatorFree :: (Coord t)
               => NArray i t          -- ^ input array
               -> [Name]              -- ^ row dimensions
               -> (Matrix t, [Name])  -- ^ (result, column dimensions)
matrixatorFree t nr = (reshape s (coords q), nc) where
    q = tridx nr t
    nc = drop (length nr) (map iName (dims q))
    s = product (map (flip size t) nc)

-- | Create a list of the substructures at the given level.
parts :: (Coord t) 
      => NArray i t
      -> Name        -- ^ index to expand
      -> [NArray i t]
parts a name | name `elem` (namesR a) = map (reorder orig) (partsRaw a name)
             | otherwise = error $ "parts: " ++ show name ++ " is not a dimension of "++(show $ namesR a)
    where orig = namesR a \\ [name]

partsRaw a name = map f (toRows m)
    where (_:ds,m) = firstIdx name a
          f t = A {dims=ds, coords=t}

tridx [] t = t
tridx (name:rest) t = A (d:ds) (join ts) where
    d = case lastIdx name t of
            ((_,d':_),_) -> d'
            _ -> error "wrong index sequence to reorder"
    ps = map (tridx rest) (partsRaw t name)
    ts = map coords ps
    ds = dims (head ps)

-- | Change the internal layout of coordinates.
-- The array, considered as an abstract object, does not change.
reorder :: (Coord t) => [Name] -> NArray i t -> NArray i t
reorder ns b | ns == namesR b = b
             | sort ns == sort (namesR b) = tridx ns b
             | otherwise = error $ "wrong index sequence " ++ show ns
                                    ++ " to reorder "++(show $ namesR b)

----------------------------------------------------------------------

-- | Apply a function (defined on hmatrix 'Vector's) to all elements of a structure.
-- Use @mapArray (mapVector f)@ for general functions.
mapArray :: Coord b => (Vector a -> Vector b) -> NArray i a -> NArray i b
mapArray f t
    | null (dims t) = scalar (f (coords t)@>0)
    | otherwise = mkNArray (dims t) (f (coords t))

liftNA2 f (A d1 v1) (A _d2 v2) = A d1 (f v1 v2)

-- | Class of compatible indices for contractions.
class (Eq a, Show (Idx a)) => Compat a where
    compat :: Idx a -> Idx a -> Bool
    opos   :: Idx a -> Idx a



contract1 t name1 name2 | ok = foldl1' (liftNA2 (+)) y
                        | otherwise = error $ "wrong contraction1: "
                                    ++(show $ dims t)++" "
                                    ++ name1++" "++name2
    where ok = (compat <$> getName t name1 <*> getName t name2) == Just True
          x = map (flip partsRaw name2) (partsRaw t name1)
          y = map head $ zipWith drop [0..] x

getName t name = d where
    l = filter ((==name).iName) (dims t)
    d = if null l
            then Nothing
            else Just (head l)

contract1c t n = contract1 renamed n n'
    where n' = " "++n++" " -- forbid spaces in names...
          renamed = renameSuperRaw (t) auxnames
          auxnames = h ++ (n':r)
          (h,_:r) = break (==n) (namesR t)

common1 t = [ n1 | (a,n1) <- x , (b,n2) <- x, a>b, n1==n2]
    where x = zip [0 ::Int ..] (namesR t)

contract t = foldl' contract1c t (common1 t)

-------------------------------------------------------------

-- | Check if two arrays have the same structure.
sameStructure :: (Eq i) => NArray i t1 -> NArray i t2 -> Bool
sameStructure a b = sortBy (compare `on` iName) (dims a) == sortBy (compare `on` iName) (dims b)

-------------------------------------------------------------

-- | Apply an element-by-element binary function to the coordinates of two arrays. The arguments are automatically made conformant.
zipArray :: (Coord a, Coord b, Compat i)
   => (Vector a -> Vector b -> Vector c) -- ^ transformation
   -> NArray i a
   -> NArray i b
   -> NArray i c
zipArray o a b = liftNA2 o a' b' where
    (a',b') = makeConformantT (a,b)

-------------------------------------------------------

-- | Create an array from a list of subarrays. (The inverse of 'parts'.)
newIndex:: (Coord t, Compat i) =>
     i  -- ^ index type
     -> Name
     -> [NArray i t]
     -> NArray i t
newIndex i name ts = r where
    ds = Idx i (length ts) name : (dims (head cts))
    cts = makeConformant ts
    r = mkNArray ds (join $ map coords cts)

-------------------------------------------------------

-- | Obtain a canonical base for the array.
basisOf :: Coord t => NArray i t -> [NArray i t]
basisOf t = map (dims t `mkNArray`) $ toRows (ident . dim . coords $ t)

-------------------------------------------------------------

instance (Coord t, Coord (Complex t), Compat i, Container Vector t) => Container (NArray i) t where
    toComplex (r,c) = zipArray (curry toComplex) r c
    fromComplex t = let (r,c) = fromComplex (coords t)
                     in (mapArray (const r) t, mapArray (const c) t)
    comp = mapArray comp
    conj = mapArray conj
    real = mapArray real
    complex = mapArray complex


-- instance (NFData t, Element t) => NFData (NArray i t) where
--     rnf = rnf . coords

----------------------------------------------------------------------

-- | obtains the common value of a property of a list
common :: (Eq a) => (b->a) -> [b] -> Maybe a
common f = commonval . map f where
    commonval :: (Eq a) => [a] -> Maybe a
    commonval [] = Nothing
    commonval [a] = Just a
    commonval (a:b:xs) = if a==b then commonval (b:xs) else Nothing

------------------------------------------------------------------------

-- | Extract the 'Matrix' corresponding to a two-dimensional array,
-- in the rows,cols order.
asMatrix :: (Coord t) => NArray i t -> Matrix t
asMatrix a | order a == 2 = reshape c (coords a')
           | otherwise = error $ "asMatrix requires a 2nd order array."
    where c = size (last (namesR a')) a'
          a' = reorder (sort (namesR a)) a

-- | Extract the 'Vector' corresponding to a one-dimensional array.
asVector :: (Coord t) => NArray i t -> Vector t
asVector a | order a == 1 = coords a
           | otherwise = error $ "asVector requires a 1st order array."

-- | Extract the scalar element corresponding to a 0-dimensional array.
asScalar :: (Coord t) => NArray i t -> t
asScalar a | order a == 0 = coords a @>0
           | otherwise = error $ "asScalar requires a 0th order array."

------------------------------------------------------------------------

-- | Create a 1st order array from a 'Vector'.
fromVector :: Compat i => i -> Vector t -> NArray i t
fromVector i v = mkNArray [Idx i (dim v) "1"] v

-- | Create a 2nd order array from a 'Matrix'.
fromMatrix :: (Compat i, Coord t) => i -> i -> Matrix t -> NArray i t
fromMatrix ir ic m = mkNArray [Idx ir (rows m) "1",
                               Idx ic (cols m) "2"] (flatten m)

------------------------------------------------------------------------

-- | Select some parts of an array, taking into account position and value.
extract :: (Compat i, Coord t)
        => (Int -> NArray i t -> Bool)
        -> Name
        -> NArray i t
        -> NArray i t
extract f name arr = reorder (namesR arr)
                   . newIndex (typeOf name arr) name
                   . map snd . filter (uncurry f)
                   $ zip [1..] (parts arr name)

-- | Apply a list function to the parts of an array at a given index.
onIndex :: (Coord a, Coord b, Compat i) =>
     ([NArray i a] -> [NArray i b])
     -> Name
     -> NArray i a
     -> NArray i b
onIndex f name t = r where
     r = if sort (namesR x) == sort (namesR t)
            then reorder (namesR t) x
            else x
     x = newIndex (typeOf name t) name (f (parts t name))

------------------------------------------------------------------------

extend alldims (A d v) = reorder (allnames) s where
    allnames = map iName alldims
    pref = alldims \\ d
    n = product (map iDim pref)
    s = A (pref++d) (join (replicate n v))

-- | Obtains most general structure of a list of dimension specifications
conformable :: Compat i => [[Idx i]] -> Maybe [Idx i]
conformable ds | ok        = Just alldims
               | otherwise = Nothing
    where alldims = nub (concat ds)
          allnames = map iName alldims
          ok = length (allnames) == length (nub allnames)

-- | Converts a list of arrays to a common structure.
makeConformant :: (Coord t, Compat i) => [NArray i t] -> [NArray i t]
makeConformant ts =
    case conformable (map dims ts) of
        Just alldims -> map (extend alldims) ts
        Nothing -> error $ "makeConformant with inconsistent dimensions "
                         ++ show (map dims ts)

-- the same version for tuples with possibly different element types
makeConformantT (t1,t2) =
    case conformable [dims t1, dims t2] of
        Just alldims -> (extend alldims t1, extend alldims t2)
        Nothing -> error $ "makeConformantT with inconsistent dimensions "
                         ++ show (dims t1, dims t2)

---------------------------------------------

takeDiagT :: (Compat i, Coord t) => NArray i t -> [t]
takeDiagT t = map (asScalar . atT t) cds where
    n = minimum (sizesR t)
    o = order t
    cds = map (replicate o) [0..n-1]

atT :: (Compat i, Coord t) => NArray i t -> [Int] -> NArray i t
atT t c = atT' c t where
    atT' cs = foldl1' (.) (map fpart cs)
    fpart k q = parts q (head (namesR q)) !! k

----------------------------------------------

-- not very smart...

-- | This is equivalent to the regular 'product', but in the order that minimizes the size of the
-- intermediate factors.
smartProduct :: (Coord t, Compat i, Num (NArray i t)) => [NArray i t] -> NArray i t
smartProduct [] = 1
smartProduct [a] = a
smartProduct [a,b] = a*b
smartProduct ts = r where
    n = length ts
    ks = [0 .. n-1]
    xs = zip ks ts
    g a b = case analyzeProduct a b of
              Nothing -> error $ "inconsistent dimensions in smartProduct: "++(show $ dims a)++" and "++(show $ dims b)
              Just (_,c) -> c
    pairs = [ ((i,j), g a b) | (i,a) <- init xs, (j,b) <- drop (i+1) xs ]
    (p,q) = fst $ minimumBy (compare `on` snd) pairs
    r = smartProduct (ts!!p * ts!!q : (dropElemPos p . dropElemPos q) ts)

dropElemPos k xs = take k xs ++ drop (k+1) xs

----------------------------------------------

-- | sequence of n indices with given prefix
seqIdx :: Int -> String -> [Name]
seqIdx n prefix = [prefix ++ show k | k <- [1 .. n] ]
