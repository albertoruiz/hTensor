import Numeric.LinearAlgebra.Exterior
import Numeric.LinearAlgebra.Array.Util(formatFixed,asMatrix)
import Numeric.LinearAlgebra (det,(><))

printAS = print . asMultivector

sh = putStrLn . formatFixed 2

-- 'listTensor' specialized for Tensor Double
infixl 9 #
(#) :: [Int] -> [Double] -> Tensor Double
(#) = listTensor

m = [3,-3]#[ 1,2,5,
             1,2,8,
            -2,0,4]

eps = leviCivita 3

a = vector [1,0,0] /\ vector [0,1,0]
b = vector [2,0,100,0] /\ vector [0,3,0,0] /\ vector [0,0,4,0]

im = eps!"ijb"* m!"pi" * m!"qj" * cov eps!"apq"

main = do
    putStrLn "exterior product"
    print a
    sh a
    printAS a
    printAS b
    putStrLn "\ndeterminant"
    printAS $ cov eps!"pqr" * m!"pi" * m!"qj" * m!"rk"
    print $ det (asMatrix m)
    putStrLn "\ninverse"
    sh $ im
    sh $ im!"ik" * m!"kj"
    putStrLn "\nmeet and join"
    printAS $ (vector [1,0,1] /\ vector [0,1,0]) \/ (vector [1,1,0] /\ vector [0,0,1])
    putStrLn "\nEuclidean inner product of r-vectors"
    printAS $ (vector [1,0,1] /\ vector [0,1,0]) `inner` (vector[1,1,1])
    print $ vector[3,5] `inner` vector[2,1]
