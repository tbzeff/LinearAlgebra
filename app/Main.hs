module Main (main) where

import Matrices
import LAAlgorithms

main :: IO ()
main = testQR

testQR :: IO ()
testQR = do 
    let mat :: (RealFrac a, Floating a) => Mat a
        mat = toMat [[1,1,0],[1,0,1],[0,1,1]]
    let qr :: (RealFrac a, Floating a) => (Mat a, Mat a)
        qr = qrDecomp mat
    let q :: (RealFrac a, Floating a) => Mat a
        q = fst qr
    let r :: (RealFrac a, Floating a) => Mat a
        r = snd qr
    putStrLn ""
    putStrLn $ "Initial Matrix A ="
    prettyPrintMat mat
    putStrLn "Matrix Q ="
    prettyPrintMat q 
    putStrLn "Matrix R ="
    prettyPrintMat r
    let multiplied :: (RealFrac a, Floating a) => Mat a
        multiplied = applyPrecisionMat (10**(-12)) $ mmMult r q
    if multiplied == mat
    then do 
        putStrLn "Q*R = A"
        putStrLn $ "Determinant: " ++ (show $ determinant mat)
    else do 
        putStrLn "Q*R /= A"
        putStrLn $ "Q*R = " ++ (show multiplied)
    putStrLn ""
    putStrLn "Eigenvalues from QR Algorithm:"
    prettyPrintMat $ qrAlgorithm mat
        