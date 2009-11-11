import Numeric.LinearAlgebra.Array
import Numeric.LinearAlgebra.Array.Util
import Numeric.LinearAlgebra.Array.Decomposition
import Numeric.LinearAlgebra hiding ((.*))
import System.Random
import Text.Printf

-- 'listArray' specialized for Array Double
infixl 9 #
(#) :: [Int] -> [Double] -> Array Double
(#) = listArray

sh x = putStrLn . formatFixed 2 $ x

-----------------------------------------------

hard = [10,2,3,4]#[1..]


randomArray dims = do
    seed <- randomIO
    let g = mkStdGen seed
        cs = randomRs (-1,1) g
    return (dims # cs)


randomArrayRank r ns = do
    seed <- randomIO
    let dummy = ns # [0,0..]
        d0:as0 = cpInitRandom seed dummy r
        as = map f as0 where f x = (map unitT `onIndex` (head (names x))) x
        d = d0 .* diagT (map fromIntegral [r, r-1 .. 1]) (length ns) `rename` (names d0)
    return (d:as)

-----------------------------------------------

separ msg = putStrLn (replicate 60 '=') >> putStrLn msg

---------------------------------------------------------------

cp finit delta epsilon t = cpAuto (finit t) delta epsilon t
iS   = cpInitSvd . fst . hosvd
iR s = cpInitRandom s

analyzeCP finit t delta epsilon = do
    putStrLn $ "\nCP delta=" ++ show delta ++ ", epsilon=" ++ show epsilon
    putStr "dims = "; print (sizes t)
    let y = cp finit delta epsilon t
    print (takeDiagT (head y))
    printf "error = %.3f %% \n" $ 100 * frobT (t - product y) / frobT t

frobT = pnorm PNorm2 . coords

unitT t = t / scalar (frobT t)


testCP finit r ns delta epsilon = do
    t <- product `fmap` randomArrayRank r ns
    analyzeCP finit t delta epsilon

main = do
    separ "rank 2, random init"
    testCP (iR 1000) 2 [3,4,5] 0.01 1E-6

    separ "rank 4, random init"
    testCP (iR 500) 4 [2,4,5] 0.01 1E-6

    separ "rank 5, hosvd init"
    testCP iS 5 [10,11,12,13,14] 0.01 1E-6

    separ "rank 3, fixed rank 3 approximation"
    t <- product `fmap` randomArrayRank 3 [3,4,5]
    let (z,errs) = cpRun (cpInitSvd (fst $ hosvd t) 3) 0.01 1E-6 t
    print $ takeDiagT (head z)
    print (head errs)

    separ "rank 3, fixed rank 2 approximation"
    t <- product `fmap` randomArrayRank 3 [3,4,5]
    let (z,errs) = cpRun (cpInitSvd (fst $ hosvd t) 2) 0.01 1E-6 t
    print $ takeDiagT (head z)
    print (head errs)

    separ "'increasing' coordinates, hosvd init"
    let z = cp iS 0.01 1E-6 hard
    print $ takeDiagT (head z)
    sh (product z)

    separ "random coordinates, random init"
    t <- randomArray [3,4,5]
    analyzeCP (iR 300) t 1 1

    separ "random coordinates, random init"
    t <- randomArray [3,4,5,6]
    analyzeCP (iR 300) t 10 1
