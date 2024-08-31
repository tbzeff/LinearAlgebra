module Matrices (
    Vec (Vec), Mat (Mat), toVec, toMat, mapVec, numCols, numRows, isZero, rotateUp, 
    rotateDown, swapRows, scaleRow, scale, applyPrecisionThreshold, 
    applyPrecisionVec, applyPrecisionMat, magnitude, normalize, vSubt, vAdd,
    getCol, getRow, vgetElem, mgetElem, mvMult, mmMult, transposeMat, unwrapVec,
    unwrapMat, norm1, norm2, normi, dotp, corr, proj, deleteCol, deleteRow,
    prettyPrintVec, prettyPrintMat, cofactor, minor, determinant
) where

import Data.Char (intToDigit)
import Data.List (transpose)

data Vec a = Vec [a] deriving (Eq, Show)
data Mat a = Mat [Vec a] deriving (Eq, Show)

-- Function to pretty print a Vec a
prettyPrintVec :: Show a => Vec a -> String
prettyPrintVec (Vec xs) = unwords (map show xs)

-- Function to pretty print a Mat a
prettyPrintMat :: Show a => Mat a -> IO ()
prettyPrintMat m = putStrLn $ unlines (map prettyPrintVec tm)
    where (Mat tm) = transposeMat m

toVec :: [a] -> Vec a
toVec v = Vec v

toMat :: [[a]] -> Mat a
toMat m = Mat (map toVec m)

-- Helper function to map over a Vec
mapVec :: (a -> b) -> Vec a -> Vec b
mapVec f (Vec u) = Vec (map f u)

numCols :: Mat a -> Int
numCols (Mat m) = length m

numRows :: Mat a -> Int
numRows (Mat m) = length (case m of (Vec r:_) -> r; _ -> [])

isZero :: (Num a, Ord a) => Vec a -> Bool
isZero (Vec []) = True
isZero (Vec (x:xs))
    | 0 < abs x = False
    | otherwise = isZero (Vec xs)

rotateUp :: Vec a -> Vec a
rotateUp (Vec []) = Vec []
rotateUp (Vec (x:xs)) = Vec (xs ++ [x])

rotateDown :: Vec a -> Vec a
rotateDown (Vec []) = Vec []
rotateDown (Vec v) = Vec ((last v) : (take ((length v) - 1) v))

swapRows :: Num a => Int -> Int -> Mat a -> Mat a
swapRows i j m
    | i == j = m
    | i < 0 || (length tm) <= j = m 
    | otherwise =  transposeMat swapped
    where (Mat tm) = transposeMat m
          front = take i tm
          mid = take (j - i) (drop i tm)
          end = drop j tm
          swapped = Mat (front ++ [getRow j m] ++ mid ++ [getRow i m] ++ end)

scaleRow :: Fractional a => Int -> a -> Mat a -> Mat a
scaleRow i c m = if i < 1 || (length tm) < i
                 then m
                 else transposeMat scaled
    where (Mat tm) = transposeMat m
          scaled = Mat ((take (i - 1) tm) ++ [scale c (getRow (i - 1) m)] ++ (drop i tm))

-- Function to scale a vector
scale :: Fractional a => a -> Vec a -> Vec a
scale c u = mapVec (*c) u

applyPrecisionThreshold :: (RealFrac a) => a -> a -> a
applyPrecisionThreshold threshold x 
    | roundedCondition = rounded
    | otherwise = x
    where rounded = fromIntegral $ round x
          roundedCondition = 
            if (abs x) - (abs rounded) < 0 
            then (abs rounded) - threshold < (abs x) 
            else (abs rounded) + threshold > (abs x)

applyPrecisionVec :: (RealFrac a) => a -> Vec a -> Vec a
applyPrecisionVec threshold v = mapVec (\x -> applyPrecisionThreshold threshold x) v

applyPrecisionMat :: (RealFrac a) => a -> Mat a -> Mat a
applyPrecisionMat threshold (Mat m) = Mat $ map (\v -> applyPrecisionVec threshold v) m

magnitude :: (Fractional a, Floating a) => Vec a -> a
magnitude (Vec v) = sqrt . sum $ map (^2) v

normalize :: (Fractional a, Floating a) => Vec a -> Vec a
normalize v = mapVec (*(1/(magnitude v))) v

vSubt :: Num a => Vec a -> Vec a -> Vec a
vSubt (Vec m) (Vec n) = Vec (zipWith (-) m n)

vAdd :: Num a => Vec a -> Vec a -> Vec a
vAdd (Vec m) (Vec n) = Vec (zipWith (+) m n)

-- Get column i from a matrix
getCol :: Int -> Mat a -> Vec a
getCol i (Mat a) = (a !! i)

-- Get row j from a matrix
getRow :: Int -> Mat a -> Vec a
getRow j m = getCol j $ transposeMat m

deleteCol :: Eq a => Int -> Mat a -> Mat a
deleteCol i (Mat m) = Mat $ filter (/= (getCol i (Mat m))) m

deleteRow :: Eq a => Int -> Mat a -> Mat a
deleteRow j (Mat m) = transposeMat $ Mat $ filter (/= (getRow j (Mat m))) tm
    where (Mat tm) = transposeMat (Mat m)

-- Get element at index i from a vector
vgetElem :: Int -> Vec a -> a
vgetElem i (Vec u) = u !! i

-- Get element at (i, j) from a matrix
mgetElem :: Int -> Int -> Mat a -> a
mgetElem i j mat = vgetElem i (getCol j mat)

-- Matrix-vector multiplication
mvMult :: (RealFrac a, Floating a) => Vec a -> Mat a -> Vec a
mvMult vec m = Vec [dotp vec row | row <- tm]
    where (Mat tm) = transposeMat m

-- Matrix-matrix multiplication
mmMult :: (RealFrac a, Floating a) => Mat a -> Mat a -> Mat a
mmMult (Mat a) matB = Mat [mvMult col matB | col <- a]

-- Vector-matrix product (recursive implementation)
_vmprod :: (RealFrac a, Floating a) => Vec a -> Mat a -> [a]
_vmprod (Vec a) (Mat []) = []
_vmprod (Vec a) (Mat (x:xs)) = (dotp (Vec a) x) : (_vmprod (Vec a) (Mat xs))

-- Vector-matrix product (wrapper function)
vmprod_ :: (RealFrac a, Floating a) => Vec a -> Mat a -> [a]
vmprod_ vec mat = _vmprod vec mat

-- Transpose a matrix
transposeMat :: Mat a -> Mat a
transposeMat (Mat a) = Mat (map Vec (transpose (map unwrapVec a)))

-- Helper function to unwrap a vector
unwrapVec :: Vec a -> [a]
unwrapVec (Vec a) = a

unwrapMat :: Mat a -> [[a]]
unwrapMat (Mat a) = map unwrapVec a

linearlyIndependent :: Mat a -> Bool
linearlyIndependent = undefined

-- Compute L1 norm of a vector
norm1 :: Num a => Vec a -> a
norm1 (Vec u) = sum (map abs u)

-- Compute L2 norm of a vector
norm2 :: (RealFrac a, Floating a) => Vec a -> a
norm2 (Vec u) = sqrt (sum (map (^2) u))

-- Compute infinity norm of a vector
normi :: (Ord a, Num a) => Vec a -> a
normi (Vec u) = maximum (map abs u)

-- Dot product of two vectors
dotp :: (RealFrac a, Floating a) => Vec a -> Vec a -> a
dotp (Vec u) (Vec v) = sum $ zipWith (*) u v

-- Correlation between two vectors
corr :: (RealFrac a, Floating a) => Vec a -> Vec a -> a
corr vec1 vec2 = (dotp vec1 vec2) / ((norm2 vec1) * (norm2 vec2))

-- Projection of vector v onto vector u
proj :: (RealFrac a, Floating a) => Vec a -> Vec a -> Vec a
proj vecV vecU = scale scalar vecU
    where scalar = (dotp vecV vecU) / (dotp vecU vecU)

minor :: (Floating a, Eq a) => Int -> Int -> Mat a -> a
minor i j m = determinant $ deleteRow j $ deleteCol i m

cofactor :: (Floating a, Eq a, Floating a) => Int -> Int -> Mat a -> a
cofactor i j m = ((-1) ** fromIntegral (i + j)) * (minor i j m)

det2x2 :: (Num a) => Mat a -> a
det2x2 m = ((mgetElem 0 0 m) * (mgetElem 1 1 m)) - ((mgetElem 1 0 m) * (mgetElem 0 1 m))

determinant :: (Floating a, Eq a) => Mat a -> a
determinant (Mat [Vec [x]]) = x  -- Base case for 1x1 matrix
determinant mat@(Mat m)
    | numCols mat == 2 = det2x2 mat
    | otherwise = sum [cofactor i 0 mat * mgetElem 0 i mat | i <- [0..numCols mat - 1]]

-- increase the "cell depth" of a matrix
increaseDepth :: (Eq a) => Mat a -> Mat a
increaseDepth (Mat m) = transposeMat finish
    where 
        expCols [] = []
        expCols (x:xs) = x : x : expCols xs
        (Mat expRows) = transposeMat (Mat (expCols m))
        finish = Mat $ expCols expRows

decreaseDepth :: (Eq a) => Mat a -> Mat a
decreaseDepth (Mat m) = transposeMat reducedCols
    where
        reduceCols [] = []
        reduceCols (x:y:xs) = x : reduceCols xs
        (Mat reducedRows) = transposeMat (Mat (reduceCols m))
        reducedCols = Mat $ reduceCols reducedRows