-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Packed.Array.Display
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  provisional
-- Portability :  portable
--
-- Formatting utilities

-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Array.Display (
    formatArray, formatFixed, formatScaled, printA, dummyAt, noIdx, showBases,
) where

import Numeric.LinearAlgebra.Array.Internal
import Data.Packed
import Data.List
import Text.Printf

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
        | order t == 0 = rect (f (coords t @> 0))
        | order t == 1 =  rect $ unwords $ map f (toList $ coords t)
        | order t == 2 =  decor t $ rect $ w1 $ format " " f (reshape (iDim $ last $ dims t) (coords t))
        | otherwise    = decor t (g ps)
      where ps = map (fmt gs ) (partsRaw t (head (names t)))
    ds = showNice (filter ((/='*').head.iName) $ dims x)
    addds = if null ds then (showRawDims (dims x) :) else (ds:)
    w1 = unlines . map (' ':) . lines
    ms = cycle [dispV 1, dispH 2]
    decor t | odd (order t) = id
            | otherwise = decorLeft  (names t!!0) . decorUp (names t!!1)


showNice x = unwords . intersperse "x" . map show $ x
showRawDims = showNice . map iDim . filter ((/="*").iName)

------------------------------------------------------

-- | Show a multidimensional array as a nested 2D table.
formatArray :: (Coord t, Compat i)
      => (t -> String) -- ^ format function (eg. printf \"5.2f\")
      -> NArray i t
      -> String
formatArray f t | odd (order t) = formatAux f (dummyAt 0 t)
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

-- | Insert a dummy index of dimension 1 at a given level (for formatting purposes).
dummyAt :: Int -> NArray i t -> NArray i t
dummyAt k t = mkNArray d' (coords t) where
    (d1,d2) = splitAt k (dims t)
    d' = d1 ++ d : d2
    d = Idx 1 "*" (iType (head (dims t))) -- undefined

-- | Rename indices so that they are not shown in formatted output.
noIdx :: Compat i => NArray i t -> NArray i t
noIdx t = renameRaw t (map ('*':) (names t))
