import Numeric.LinearAlgebra.Array
import Numeric.LinearAlgebra.Array.Util
import Numeric.LinearAlgebra.Array.Solve
import Numeric.LinearAlgebra.Array.Decomposition
import Numeric.LinearAlgebra hiding ((.*))
import Data.List
import System.Random
import Text.Printf
import Control.Monad

infixl 9 #
(#) :: [Int] -> [Double] -> Array Double
(#) = listArray

sh x = putStrLn . formatFixed 2 $ x

mk ds f = ds # map f (sequence $ map (enumFromTo 1 . fromIntegral) $ ds)

withOnesAt ds ps = mk ds f where
    f p = if p `elem` ps then 1 else 0

randomArray dims = do
    seed <- randomIO
    let g = mkStdGen seed
        cs = randomRs (-1,1) g
    return (dims # cs)

----------------------------------------------

core = [4,13,4] `withOnesAt` [
    [1,1,1],
    [1,2,2],
    [1,3,3],
    [1,10,4],
    [2,4,1],
    [2,5,2],
    [2,6,3],
    [2,11,4],
    [3,7,1],
    [3,8,2],
    [3,9,3],
    [3,12,4],
    [4,13,4]]!"CFX"

------------------------------------------------

main = do
    frames <- randomArray [13,30]
    points <- randomArray [4,20]
    cameras <- randomArray [4,15]
    let obs = core * cameras!"Cc" * frames!"Ff" * points!"Xx" 
    sh $ dummyAt 1 obs
    let (f,ss) = hosvd obs
--    mapM_ print ss
    mapM_ (print . rank . flip fibers obs ) (names obs)
    let q = truncateFactors [4,13,4] f
{-    sh (dummyAt 1 $ head f)
    sh (dummyAt 1 $ head q)
    sh (dummyAt 1 $ core)-}
--     mapM_ sh $ als 1 1 core (head q)
--     print (dims core)
--     print $ dims (product $ als 1 1 core (head q))
    let (sol,es) = multilinearSolve eqnorm 0.1 1E-3 (head q) core
    print (head es)
    print (length es)
    print $ length $ snd $ multilinearSolve id 0.1 1E-3 (head q) core
 --   let s0 = als0 core (head q)
 --   mapM_ sh s0
 --   let sol = alsArg core 3 $ alsArg core 3 $ alsArg core 2 $ alsArg core 1 s0

 --   mapM_ sh sol
    let s = concat [sol,tail q]
    mapM_ (print.dims) s
    print (100 * frobT (product s - obs)/frobT obs)
    sh (head s)
    let cam = s!!4 * s!!1
        frm = s!!5 * s!!2
        pts = s!!6 * s!!3
    putStrLn "----------------------"
--    sh cam
--    sh frm
    sh pts
    print (100 * frobT (product [core,cam,frm,pts] - obs)/frobT obs)
    print $ map frobT sol

-- infixl 5 ./
-- (./) = zipArray (/)

frobT = pnorm PNorm2 . coords
unitT t = t * scalar (recip (frobT t))

-----------------------------

-- a = [2,4,3,2]#[1..] !"pqrs"
-- b = [3,5,7,4]#[5,9..] !"ijkl"

b = [2,4,3,2]#[1..] !"pqrs"
a = [3,5,4,4]#[5,9..] !"ijkl"


--a = [2,2,2,2]#[1..] !"pqrs"
--b = [3,3,3,3]#[5,9..] !"ijkl"


kk = multilinearSolve eqnorm 0.1 1E-3  a b

s0 = alsInitRandom 170 a b

pru f a b = head $ snd $ multilinearSolve f 0.1 1E-3  a b

test = do
    sh a
    sh b
    print (head $ snd kk)
    print (length (snd kk))
    sh (product (fst kk))
    mapM_ sh (fst kk)
    print $ map frobT (fst kk)


p0 = alsArg a 2 $ alsArg a 1 s0
q0 = alsArg' a 2 $ alsArg' a 1 s0
p = alsArg a 3 p0
q = alsArg' a 3 p0
--p = alsArg a 3 s0
--q = alsArg' a 3 s0
xx = alsArg' a 1 $ alsArg' a 4 $ alsArg' a 3 $ alsArg' a 2 $ alsArg' a 1 s0

test2 = do
    
--    sh $ unitT $ p0!!3
--    sh $ unitT $ q0!!3
    sh $ id $ p!!3
    sh $ id $ q!!3
    print $ 100 * frobT (p!!3 - q!!3) / frobT (p!!3)
--    mapM_ sh $ tail $ alsArg a 3 s0
--    mapM_ sh $ tail $ alsArg' a 3 s0
--    print $ frobT (product p0 - product q0)
--    mapM_ (print . dims) s0
--    mapM_ (print . dims) p0
--    mapM_ (print . dims) q0
--    mapM_ (print . dims) p
--    mapM_ (print . dims) q
--    mapM_ (print . dims) xx
