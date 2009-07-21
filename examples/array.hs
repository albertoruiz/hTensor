import Numeric.LinearAlgebra.Array
import Data.Packed.Array.Util
import Control.Applicative
import Data.List
import Text.Printf
import Data.Packed

-- 'listArray' specialized for Array Double
infixl 9 #
(#) :: [Int] -> [Double] -> Array Double
(#) = listArray

(<|) :: Name -> [Array Double] -> Array Double
infixl 8 <|
n <| ls = index n ls

i = ("i" <|)
j = ("j" <|)
k = ("k" <|)

sh x = printS 2 x

a = [3,4,2] # [1..]
b = [2,3] # [5,6]
c = 7 :: Array Double
s = [2,2,2,2] # [1..]

t = [3,3,3,3]#[1 ..]!"ijkl"

q = [2,4,3] # (fun <$> r 2 <*> r 4 <*> r 3) !"ijk"
    where r k = [1..k]
          fun = \i j k -> i*2*j-k

m = j [i[2,0,0], i[1,0,1], i[0,3,0]] ~> "ij"

main = do
    putStrLn "8-dimensional array"
    sh $ (replicate 8 2)#[1::Double ..]!(take 8 ['a'..])
    ------------------------
    putStrLn "different display formats"
    sh $ a!"ijk"
    printA "%7.3f" a
    putStrLn . formatS 2 $ a!"ijk"
    ------------------------
    putStrLn "array defined using a function"
    sh q
    ------------------------
    putStrLn "contraction"
    sh t
    sh $ t!"ijkk"
    ------------------------
    putStrLn "tensor product"
    sh $ m
    sh $ (t !"pqrs" * m!"kr") ~> "pqks"
    ------------------------
    putStrLn "automatic conformability"
    sh $ j[1,2,3] + k[10,20]
    sh $ k [m, 3*m-1, 7, i[1,2,3]]
