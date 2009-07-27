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
    rank, names, size, typeOf , dims, coords,
    Compat(..),
    -- * Array creation
    scalar,
    mkNArray,
    fromVector, fromMatrix,
    -- * Array manipulation
    rename,(!),
    parts,
    (|*|),
    zipArray,
    mapArray,
    extract,
    onIndex,
    -- * Utilities
    reorder, (~>),
    sameStructure,
    conformable,
    makeConformant,
    mapTypes,
    renameRaw,
    formatArray, formatFixed, formatScaled, printA,
    showBases,
    newIndex,
    dummyAt, noIdx,
    basisOf,
    common,
    Coord,
    asMatrix, asVector, asScalar
) where

import Data.Packed
import Data.List
import Numeric.LinearAlgebra(outer,multiply,Field)
import Control.Applicative
import Data.Function(on)
import Text.Printf

-- | Types that can be elements of the multidimensional arrays.
class (Num (Vector t), Field t) => Coord t
instance Coord Double
instance Coord (Complex Double)

-- import Debug.Trace
-- 
-- debug s f x = trace (s ++ ": " ++ show (f x)) x

-- | indices are denoted by strings, (frequently single-letter)
type Name = String

-- | Dimension descriptor.
data Idx i = Idx { iDim  :: Int
                 , iName :: Name
                 , iType :: i
                 } deriving (Eq)



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


-- | 'rename' the indices with single-letter names. Equal indices of compatible type are contracted out.
infixl 8 !
(!) :: (Coord t, Compat i)
       => NArray i t
       -> String   -- ^ new indices
       -> NArray i t
t ! ns = rename t (map return ns)

-- | Rename indices. Equal indices are contracted out.
rename :: (Coord t, Compat i)
       => NArray i t
       -> [Name]     -- ^ new names
       -> NArray i t
rename t ns = reorder orig (contract t')
    where t' = renameRaw t ns
          orig = nub (names t') \\ common1 t'


renameRaw (A d v) l | length l == length d = A d' v
                    | otherwise = error $ "rename " ++ show d ++ " with " ++ show l
    where d' = zipWith f d l
          f i n = i {iName=n}

mapDims f (A d v) = A (map f d) v

mapTypes :: (i1 -> i) -> NArray i1 t -> NArray i t
mapTypes f = mapDims (\i -> i {iType = f (iType i)})

-- mapNames f = mapDims (\i -> i {iName = f (iName i)})

-- | Index names.
names :: NArray i t -> [Name]
names = map iName . dims

-- | Dimension of given index.
size :: Name -> NArray i t -> Int
size n t = (iDim . head) (filter ((n==).iName) (dims t))

-- | Type of given index.
typeOf :: Compat i => Name -> NArray i t -> i
typeOf n t = (iType . head) (filter ((n==).iName) (dims t))


-- | The number of dimensions of a multidimensional array.
rank :: NArray i t -> Int
rank = length . dims

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


-- | Create a list of the substructures at the given level.
parts :: (Coord t) 
      => NArray i t
      -> Name        -- ^ index to expand
      -> [NArray i t]
parts a name | name `elem` (names a) = map (reorder orig) (partsRaw a name)
             | otherwise = error $ "parts: " ++ show name ++ " is not a dimension of "++(show $ names a)
    where orig = names a \\ [name]

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
reorder ns b | ns == names b = b
             | sort ns == sort (names b) = tridx ns b
             | otherwise = error $ "wrong index sequence " ++ show ns
                                    ++ " to reorder "++(show $ names b)


-- | 'reorder' (transpose) the dimensions of the array (with single letter names).
--
-- Operations are defined by named indices, so the transposed array is operationally equivalent to the original one.
infixl 8 ~>
(~>) :: (Coord t) => NArray i t -> String -> NArray i t
t ~> ns = reorder (map return ns) t

-------------------------------------------------------------

rawProduct (A d1 v1) (A d2 v2) = A (d1++d2) (flatten (outer v1 v2))

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
          renamed = renameRaw (t) auxnames
          auxnames = h ++ (n':r)
          (h,_:r) = break (==n) (names t)

common1 t = [ n1 | (a,n1) <- x , (b,n2) <- x, a>b, n1==n2]
    where x = zip [0 ::Int ..] (names t)

contract t = foldl' contract1c t (common1 t)

----------------------------------------------------------------------

contract2 t1 t2 n | ok = A (tail ds1 ++ tail ds2) (flatten m)
                  | otherwise = error $ "wrong contraction2: "++ n ++ " of "++
                                      (show $ dims t1)++" and "++ (show $ dims t2)
  where ok = (compat <$> getName t1 n <*> getName t2 n) == Just True
        (ds1,m1) = firstIdx n t1
        (ds2,m2) = firstIdx n t2
        m = (trans m1) `multiply` m2

common2 t1 t2 = [ n1 | n1 <- names t1, n2 <- names t2, n1==n2]

infixl 5 |*|
-- | Tensor product with automatic contraction of repeated indices, following Einstein summation convention.
(|*|) :: (Coord t, Compat i)
      => NArray i t -> NArray i t -> NArray i t
t1 |*| t2 = r where
    cs = common2 t1 t2
    r = case cs of
        [] -> rawProduct t1 t2
        n:_ -> reorder orig $ contract (contract2 t1 t2 n)
    orig = nub (names t1 ++ names t2) \\ cs

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

showBases x = f $ concatMap (shbld) x
    where   shbld (c,[]) = shsign c ++ showc c
            shbld (c,l) = shsign c ++ g (showc c) ++ "{"++ concatMap show l++"}"
            shsign c = if c < 0 then " - " else " + "
            showc c
                | abs (fromIntegral (round c :: Int) - c) <1E-10  = show (round $ abs c::Int)
                | otherwise = printf "%.3f" (abs c)
            f (' ':'+':' ':rs) = rs
            f (' ':'-':' ':rs) = '-':rs
            f a = a
            g "1" = ""
            g a = a

---------------------------------------------------------

data Rect = Rect { li :: Int, co :: Int, els :: [String] }

rect s = pad r c (Rect r 0 ss)
    where ss = lines s
          r  = length ss
          c  = maximum (map length ss)

pad nr nc (Rect r c ss) = Rect (r+r') (c+c') ss'' where
    r' = max 0 (nr-r)
    c' = max 0 (nc-c)
    ss' = map (padH nc) ss
    ss'' = replicate r' (replicate nc '-') ++ ss'
    padH l s = take (l-length s) (" | "++repeat ' ') ++ s

dispH :: Int -> [Rect] -> Rect
dispH k rs = Rect nr nc nss where
    nr = maximum (map li rs)
    nss' = mapTail (\x-> pad nr (co x + k) x) rs
    nss = foldl1' (zipWith (++)) (map els nss')
    nc = length (head nss)

dispV :: Int -> [Rect] -> Rect
dispV k rs = Rect nr nc nss where
    nc = maximum (map co rs)
    nss' = mapTail (\x-> pad (li x + k) nc x) rs
    nss = concatMap els nss'
    nr = length nss

mapTail f (a:b) = a : map f b
mapTail _ x     = x




formatAux f x = unlines . addds . els . fmt ms $ x where
    fmt [] _ = undefined -- cannot happen
    fmt (g:gs) t
        | rank t == 0 = rect (f (coords t @> 0))
        | rank t == 1 =  rect $ unwords $ map f (toList $ coords t)
        | rank t == 2 =  decor t $ rect $ w1 $ format " " f (reshape (iDim $ last $ dims t) (coords t))
        | otherwise    = decor t (g ps)
      where ps = map (fmt gs ) (partsRaw t (head (names t)))
    ds = showNice (filter ((/='*').head.iName) $ dims x)
    addds = if null ds then (showRawDims (dims x) :) else (ds:)
    w1 = unlines . map (' ':) . lines
    ms = cycle [dispV 1, dispH 2]
    decor t | odd (rank t) = id
            | otherwise = decorLeft  (names t!!0) . decorUp (names t!!1)


showNice x = unwords . intersperse "x" . map show $ x
showRawDims = showNice . map iDim . filter ((/="*").iName)

------------------------------------------------------

-- | Show a multidimensional array as a nested 2D table.
formatArray :: (Coord t, Compat i)
      => (t -> String) -- ^ format function (eg. printf \"5.2f\")
      -> NArray i t
      -> String
formatArray f t | odd (rank t) = formatAux f (dummyAt 0 t)
            | otherwise    = formatAux f t


decorUp s rec
    | head s == '*' = rec
    | otherwise     = dispV 0 [rs,rec]
  where
    c = co rec
    c1 = (c - length s) `div` 2
    c2 = c - length s - c1
    rs = rect $ replicate c1 ' ' ++ s ++ replicate c2 ' '

decorLeft s rec
    | head s == '*' = rec
    | otherwise     = dispH 0 [rs,rec]
  where
    c = li rec
    r1 = (c - length s+1) `div` 2
    r2 = c - length s - r1
    rs = rect $ unlines $ replicate r1 spc ++ s : replicate (r2) spc
    spc = replicate (length s) ' '

------------------------------------------------------

-- | Print the array as a nested table with the desired format (e.g. %7.2f) (see also 'formatArray', and 'formatScaled').
printA :: (Coord t, Compat i, PrintfArg t) => String -> NArray i t -> IO ()
printA f t = putStrLn (formatArray (printf f) t)


-- | Show the array as a nested table with autoscaled entries.
formatScaled :: (Compat i)
      => Int -- ^ number of of decimal places
      -> NArray i Double
      -> String
formatScaled dec t = unlines (('(':d++")  E"++show o) : m)
    where ss = formatArray (printf fmt. g) t
          d:m = lines ss
          g x = x/10^(o::Int)
          o = floor $ maximum $ map (logBase 10 . abs) $ toList $ coords t
          fmt = '%':show (dec+3) ++ '.':show dec ++"f"

-- | Show the array as a nested table with a \"\%.nf\" format. If all entries
-- are approximate integers the array is shown without the .00.. digits.
formatFixed :: (Compat i)
      => Int -- ^ number of of decimal places
      -> NArray i Double
      -> String
formatFixed dec t
    | isInt t   = formatArray (printf ('%': show (width t) ++".0f")) t
    | otherwise = formatArray (printf ('%': show (width t+dec+1) ++"."++show dec ++"f")) t

isInt = all lookslikeInt . toList . coords
lookslikeInt x = show (round x :: Int) ++".0" == shx || "-0.0" == shx
    where shx = show x
-- needsSign t = vectorMin (coords t) < 0
-- width :: Compat i => NArray i Double -> Int
width = maximum . map (length . (printf "%.0f"::Double->String)) . toList . coords
-- width t = k + floor (logBase 10 (max 1 $ vectorMax (abs $ coords t))) :: Int
--      where k | needsSign t = 2
--              | otherwise = 1

------------------------------------------------------

-- | Create an array from a list of subarrays. (The inverse of 'parts'.)
newIndex:: (Coord t, Compat i) =>
     i  -- ^ index type
     -> Name
     -> [NArray i t]
     -> NArray i t
newIndex i name ts = r where
    ds = Idx (length ts) name i : (dims (head cts))
    cts = makeConformant ts
    r = mkNArray ds (join $ map coords cts)


-- | Insert a dummy index of dimension 1 at a given level (for formatting purposes).
dummyAt :: Int -> NArray i t -> NArray i t
dummyAt k t = mkNArray d' (coords t) where
    (d1,d2) = splitAt k (dims t)
    d' = d1 ++ d : d2
    d = Idx 1 "*" undefined

-- | Rename indices so that they are not shown in formatted output.
noIdx :: Compat i => NArray i t -> NArray i t
noIdx t = renameRaw t (map ('*':) (names t))

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
asMatrix a | rank a == 2 = reshape c (coords a)
           | otherwise = error $ "asMatrix requires a rank 2 array."
    where c = size (last (names a)) a

-- | Extract the 'Vector' corresponding to a one-dimensional array.
asVector :: (Coord t) => NArray i t -> Vector t
asVector a | rank a == 1 = coords a
           | otherwise = error $ "asVector requires a rank 1 array."

-- | Extract the scalar element corresponding to a 0-dimensional array.
asScalar :: (Coord t) => NArray i t -> t
asScalar a | rank a == 0 = coords a @>0
           | otherwise = error $ "asScalar requires a rank 0 array."

------------------------------------------------------------------------

-- | Create a rank-1 array from an hmatrix 'Vector'.
fromVector :: Compat i => i -> Vector t -> NArray i t
fromVector i v = mkNArray [Idx (dim v) "1" i ] v

-- | Create a rank-2 array from an hmatrix 'Matrix'.
fromMatrix :: (Compat i, Coord t) => i -> i -> Matrix t -> NArray i t
fromMatrix ir ic m = mkNArray [Idx (rows m) "1" ir,
                               Idx (cols m) "2" ic] (flatten m)

------------------------------------------------------------------------

-- | Select some parts of an array, taking into account position and value.
extract :: (Compat i, Coord t)
        => (Int -> NArray i t -> Bool)
        -> Name
        -> NArray i t
        -> NArray i t
extract f name arr = reorder (names arr)
                   . newIndex (typeOf name arr) name
                   . map snd . filter (uncurry f)
                   $ zip [1..] (parts arr name)

-- | Apply a list function to the parts of an array at a given index.
onIndex :: (Coord a, Coord b, Compat i) =>
     ([NArray i a] -> [NArray i b])
     -> Name
     -> NArray i a
     -> NArray i b
onIndex f name t = reorder (names t) $ newIndex (typeOf name t) name (f (parts t name))

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
