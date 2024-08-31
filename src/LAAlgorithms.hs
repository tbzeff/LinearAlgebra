module LAAlgorithms (
    orthogonalizeAndNormalize, gsProcess, qrDecomp, qrAlgorithm, gaussJordanElimination
) where

import Matrices
import Control.Monad (foldM)

-- Helper function to perform vector subtraction and normalization
orthogonalizeAndNormalize :: (RealFrac a, Floating a) => Vec a -> [Vec a] -> Vec a
orthogonalizeAndNormalize v [] = normalize v
orthogonalizeAndNormalize v (x:xs) = orthogonalizeAndNormalize (v `vSubt` (proj v x)) xs

-- Gram-Schmidt process
gsProcess :: (RealFrac a, Floating a) => Mat a -> Mat a -> Mat a
gsProcess (Mat n) (Mat []) = Mat (map normalize n)
gsProcess (Mat n) (Mat m)
    | null n = gsProcess (Mat [head m]) (Mat (tail m))
    | otherwise = gsProcess (Mat (n ++ [nextVec])) (Mat (tail m))
    where nextVec = orthogonalizeAndNormalize (head m) n

qrDecomp :: (RealFrac a, Floating a) => Mat a -> (Mat a, Mat a)
qrDecomp m = (q, r)
    where q = gsProcess (Mat []) m
          r = applyPrecisionMat (10**(-12)) $ mmMult m (transposeMat q)

-- QR algorithm to find the eigenvalues of a matrix
qrAlgorithm :: (RealFrac a, Floating a) => Int -> Mat a -> Mat a
qrAlgorithm 0 mat = mat
qrAlgorithm i mat
    | mat == newMat = mat
    | otherwise = qrAlgorithm (i - 1) newMat
    where qr = qrDecomp mat
          newMat = applyPrecisionMat (10**(-11)) $ mmMult (fst qr) (snd qr) 

qraDebug :: (RealFrac a, Floating a) => Int -> Mat a -> Mat a
qraDebug 0 mat = mat
qraDebug i mat =
    let qrd = qrDecomp mat
        q = fst qrd
        r = snd qrd

        next = applyPrecisionMat (10**(-11)) $ mmMult q r
    in qraDebug (i - 1) next


-- Recursive helper function
go :: (Fractional a, Eq a) => Mat a -> Int -> Int -> Mat a
go m r c
    | r >= (numRows m) || c >= (numCols m) = m
    | otherwise =
        case findNonZeroRow m r c of
            Nothing -> go m r (c + 1)
            Just rowIndex ->
                let m' = swapRows r rowIndex m
                    pivotRow = getRow r m'
                    pivotElement = vgetElem c pivotRow
                    m'' = if pivotElement /= 0
                          then let scaledPivotRow = scale (1 / pivotElement) pivotRow
                                   matScaled = updateRow m' r scaledPivotRow
                               in eliminateColumn matScaled r c
                          else m'
                in go m'' (r + 1) (c + 1)

-- Find the first row with a non-zero element in the current column
findNonZeroRow :: (Num a, Eq a) => Mat a -> Int -> Int -> Maybe Int
findNonZeroRow m r c
    | r >= (numRows m) = Nothing
    | mgetElem r c m /= 0 = Just r
    | otherwise = findNonZeroRow m (r + 1) c

-- Replace a row in the matrix
updateRow :: Num a => Mat a -> Int -> Vec a -> Mat a
updateRow m r newRow = transposeMat $ Mat (take r rows ++ [newRow] ++ drop (r + 1) rows)
    where (Mat rows) = transposeMat m

-- Eliminate all other entries in the current column
eliminateColumn :: Fractional a => Mat a -> Int -> Int -> Mat a
eliminateColumn m r c = foldl eliminateRow m [0..(numRows m)-1]
    where
        pivotRow = getRow r m
        eliminateRow m' i
            | i == r = m'
            | otherwise =
                let factor = mgetElem i c m'
                    newRow = vSubt (getRow i m') (scale factor pivotRow)
                in updateRow m' i newRow

-- Gauss-Jordan Elimination implementation
gaussJordanElimination :: (Fractional a, Eq a) => Mat a -> Mat a
gaussJordanElimination (Mat m) = go (Mat m) 0 0

elimColDebug :: (Fractional a, Show a) => Mat a -> Int -> Int -> IO (Mat a)
elimColDebug m r c = foldM eliminateRow m [0..(numRows m) - 1]
  where
    pivotRow = getRow r m
    eliminateRow m' i
        | i == r = return m'
        | otherwise = do
            let factor = mgetElem i c m'
            let newRow = vSubt (getRow i m') (scale factor pivotRow)
            let m'' = updateRow m' i newRow
            putStrLn $ "\nProcessing row: " ++ show i
            putStrLn $ "Pivot row (row " ++ show r ++ "): " ++ show pivotRow
            putStrLn $ "Factor (row " ++ show i ++ ", column " ++ show c ++ "): " ++ show factor
            putStrLn $ "Current row: " ++ (show (getRow i m'))
            putStrLn $ "Scaled pivot row: " ++ (show (scale factor pivotRow))
            putStrLn $ "New row " ++ show i ++ ": " ++ show newRow
            putStrLn $ "Updated matrix after eliminating row " ++ show i ++ ":"
            print m''
            return m''

goDebug :: (Fractional a, Eq a, Show a) => Mat a -> Int -> Int -> IO ()
goDebug m r c
    | r >= (numRows m) || c >= (numCols m) = do putStr "End matrix: "; print m; putStrLn ""
    | otherwise = do
        putStrLn $ "\nRow index: " ++ (show r) ++ " | Column index: " ++ (show c)
        putStr "Starting matrix: "
        print m
        case findNonZeroRow m r c of
            Nothing -> goDebug m r (c + 1)
            Just rowIndex -> do
                putStrLn $ "First non-zero row index: " ++ (show r)
                let m' = swapRows r rowIndex m
                    pivotRow = getRow r m'
                    pivotElement = vgetElem c pivotRow
                putStrLn $ "Post-Row-Swap Matrix: " ++ (show m')
                putStrLn $ "Pivot Row: " ++ (show pivotRow) ++ " | Pivot Element: " ++ (show pivotElement)
                m'' <- if pivotElement /= 0
                       then do
                            let scaledPivotRow = scale (1 / pivotElement) pivotRow
                                matScaled = updateRow m' r scaledPivotRow
                            putStrLn ("Scaled pivot row: " ++ (show scaledPivotRow))
                            putStrLn ("Scaled matrix: " ++ (show matScaled))
                            elimColDebug matScaled r c
                       else return m'
                goDebug m'' (r + 1) (c + 1)