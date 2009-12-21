-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.LinearAlgebra.Multivector
-- Copyright   :  (c) Alberto Ruiz 2009
-- License     :  GPL
--
-- Maintainer  :  Alberto Ruiz <aruiz@um.es>
-- Stability   :  experimental
--
-- A simple implementation of Geometric Algebra.
--
-- The Num instance provides the geometric product, and the Fractional
-- instance provides the inverse of multivectors.
--
-- This module provides a simple Euclidean embedding.

-----------------------------------------------------------------------------

module Numeric.LinearAlgebra.Multivector (
    Multivector, coords,
    scalar, vector, e, (/\), (-|), (\/), rever, full, rotor,
    apply,
    grade, maxGrade, maxDim,
    fromTensor
) where

import Numeric.LinearAlgebra(toList,reshape,(<\>),(@>))
import Numeric.LinearAlgebra.Array.Internal hiding (scalar,coords)
import Numeric.LinearAlgebra.Array.Display (showBases)
import Numeric.LinearAlgebra.Tensor hiding (scalar,vector)
import qualified Numeric.LinearAlgebra.Array.Internal as Array
import Data.List
import Control.Monad(filterM)
import Data.Function(on)
import qualified Data.Map as Map

powerset = filterM (const [True, False])  -- !!

base :: Int -> [[Int]]
base k = sortBy (compare `on` length) (powerset [1..k])

base' k = map (\b -> MV [(1,b)]) (base k)

data Multivector = MV { coords :: [(Double,[Int])] } deriving Eq

instance Show Multivector where
    show = showMV

maxGrade :: Multivector -> Int
maxGrade (MV l) = maximum . map (length.snd) $ l

grade :: Int -> Multivector -> Multivector
grade k (MV l) = MV $ filter ((k==).length.snd) l

maxDim :: Multivector -> Int
maxDim (MV [(_,[])]) = 0
maxDim (MV l) = maximum . concat . map snd $ l

-- | The reversion operator.
rever :: Multivector -> Multivector
rever (MV l) = MV (map r l) where
    r (c,b) = (c*fromIntegral s ,b)
        where s = signum (-1)^(k*(k-1)`div`2) :: Int
              k = length b

-- | Show the non zero coordinates of a multivector in a nicer format.
showMV :: Multivector -> String
showMV (MV x) = showBases x

-- | Creates a scalar multivector.
scalar :: Double -> Multivector
scalar s = MV [(s,[])]

-- | Creates a grade 1 multivector of from a list of coordinates.
vector :: [Double] -> Multivector
vector v = MV $ simplify $ zip v (map (:[]) [1..])


-- different product rules

-- reorders the base indices remembering the original position
r1 :: [Int] -> [(Int,[Int])]
r1 [] = []
r1 l = (m,elemIndices m l):(r1 (filter (/=m) l))
    where m = minimum l

-- geometric product
r2 :: [(Int, [Int])] -> (Double, [Int])
r2 = foldl' g (1,[])
    where g (k,l) (x,ps) = (k*s,l++t)
              where t = if even (length ps) then [] else [x]
                    s = product (map f ps')
                         where f z = if even z then 1 else -1
                               ps' = zipWith (subtract) ps [0..]

-- exterior product
r3 :: [(Int, [Int])] -> (Double, [Int])
r3 = foldl' g (1,[])
    where g (k,l) (x,ps) = (k*s,l++[x])
              where s = if length ps > 1 then 0 else if even (head ps) then 1 else -1


-- simplification and cleaning of the list of coordinates
simplify = chop . grp . sortBy (compare `on` snd)
    where grp [] = []
          grp [a] = [a]
          grp ((c1,b1):(c2,b2):rest)
              | b1 == b2  = grp ( (c1+c2,b1) : rest)
              | otherwise = (c1,b1): grp ((c2,b2):rest)
          zero (c,_) = abs c < 1E-8
          chop = cz . filter (not.zero)
          cz [] = [(0,[])]
          cz x  = x

-- sum of multivectors
gs (MV l1) (MV l2) = MV $ simplify (l1++l2)

-- geometric product
gp (MV l1) (MV l2) = MV $ simplify [g x y | x<-l1, y <-l2]
    where g (c1,b1) (c2,b2) = (k*c1*c2,b3) where (k,b3) = gpr b1 b2 --(r2.r1) (b1++b2)

-- exterior product
ge (MV l1) (MV l2) = MV $ simplify [g x y | x<-l1, y <-l2]
    where g (c1,b1) (c2,b2) = (k*c1*c2,b3) where (k,b3) = epr b1 b2 -- (r3.r1) (b1++b2)

-- contraction inner product
gi (MV l1) (MV l2) = sum [g x y | x<-l1, y <-l2]
    where g (c1,[]) (c2,is) = MV [(c1*c2,is)]
          g _ (_,[])        = 0

          g (c1,[i]) (c2,[j]) = if i==j then MV [(c1*c2,[])] else 0
          g (c1,[i]) (c2,j:js) = (g (c1,[i]) (c2,[j]) /\ MV [(1,js)])
                               - (MV [(c2,[j])] /\ g (c1,[i]) (1,js))

          g (c1,i:is) b = gi (MV [(c1,[i])]) (gi (MV[(1,is)]) (MV [b]))


instance Num Multivector where
    (+) = gs
    (*) = gp
    negate (MV l) = MV (map neg l) where neg (k,b) = (-k,b)
    abs _ = error "abs of multivector not yet defined"
    signum _ = error "signum of multivector not yet defined"
    fromInteger x = MV [(fromInteger x,[])]

instance Fractional Multivector where
    fromRational x = MV [(fromRational x,[])]
    recip (MV [(x,[])]) = MV [(recip x,[])]
    recip x = mvrecip x

-- | The k-th basis element.
e :: Int -> Multivector
e k = MV [(1,[k])]

-- | The exterior (outer) product.
(/\) :: Multivector -> Multivector -> Multivector
infixl 7 /\
(/\) = ge


-- | The contractive inner product.
(-|) :: Multivector -> Multivector -> Multivector
infixl 7 -|
(-|) = gi

-- | The full space of the given dimension. This is the leviCivita simbol, and the basis of the pseudoscalar.
full :: Int -> Multivector
full k = MV [(1,[1 .. k])] --product . map e $ [1 .. k]

-- | Intersection of subspaces.
(\/) :: Multivector -> Multivector -> Multivector
infixl 7 \/
(\/) a b = (b -| rever (full k)) -| a
    where k = max (maxDim a) (maxDim b)

-- check that it is a vector
normVec v = sqrt x where MV [(x,[])] = v * v

unitary v = v / scalar (normVec v)

-- | The rotor operator, used in a sandwich product.
rotor :: Int          -- ^ dimension of the space
      -> Double       -- ^ angle
      -> Multivector  -- ^ axis
      -> Multivector  -- ^ result
rotor k phi axis = scalar (cos (phi/2)) - scalar (sin (phi/2)) * (unitary axis*full k)


-- memoization of the rules
gprules k = Map.fromList [(x, Map.fromList [(y,(r2.r1)(x++y)) | y<-base k] )| x<-base k]

eprules k = Map.fromList [(x, Map.fromList [(y,(r3.r1)(x++y)) | y<-base k] )| x<-base k]

--reasonable limit
gpr a b = g Map.! a Map.! b
    where g = gprules 6

epr a b = g Map.! a Map.! b
    where g = eprules 6

----------------------- tensor expansion -----------------------

expand k = g
    where g (MV l) = foldl1' (zipWith (+)) $ map f l
          basepos b = m Map.! b
          m = baseraw k
          f (c,b) = en pk (basepos b) c
          pk = 2^k
          baseraw q = Map.fromList $ zip (base q) [0..]
          en n q v = replicate q 0 ++ v : replicate (n-q-1) 0

compact k t = sum $ zipWith (*) (map scalar $ toList (Array.coords t)) (base' k)

gatensor k = listTensor [-pk,-pk,pk] (concat . concat $ gacoords)
    where pk = 2^k
          gacoords = [[ f (x * y) | y<-b] | x<-b]
          b = base' k
          f = expand k

tmv k x = listTensor [2^k] (expand k x)


-- tp a b = comp (g!"ijk" * ta!"i" * tb!"j")
--     where k = max (maxDim a) (maxDim b)
--           g = gatensor k
--           ta = tmv k a
--           tb = tmv k b
--           comp = compact k

mat rowidx t = reshape c $ Array.coords t'
    where c = iDim $ last (dims t')
          t' = reorder (rowidx: (namesR t\\[rowidx])) t

-- on the right
pmat k b = mat "k" $ g!"ijk" * tb!"j"
    where g = gatensor k
          tb = tmv k b

divi k a b = compact k $ listTensor [2^k] (toList $ pmat k b <\> Array.coords (tmv k a))

mvrecip b = divi (maxDim b) 1 b

--------------------------------------------------------

-- | Extract a multivector representation from a full antisymmetric tensor.
--
-- (We do not check that the tensor is actually antisymmetric.)
fromTensor :: Tensor Double -> Multivector
fromTensor t = MV $ filter ((/=0.0).fst) $ zip vals basis
    where vals = map ((@> 0). Array.coords .foldl' partF t) (map (map pred) basis)
          r = length (dims t)
          n = iDim . head . dims $ t
          partF s i = part s (name,i) where name = iName . head . dims $ s
          basis = filter (\x-> (x==nub x && x==sort x)) $ sequence $ replicate r [1..n]

part t (name,i) = parts t name !! i

--------------------------------------------------------

-- | Apply a linear transformation, expressed as the image of the element i-th of the basis.
--
--  (This is a monadic bind!)
apply :: (Int -> Multivector) -> Multivector -> Multivector
apply f t = sum $ map g (coords t) where
    g (x,[]) = scalar x
    g (x,es) = scalar x * foldl1' (/\) (map f es)
