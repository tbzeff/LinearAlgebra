module Main (main) where

import Matrices
import LAAlgorithms

main :: IO ()
main = testQR

testQR :: IO ()
testQR = do 
    let mat :: (Integral a, RealFrac a, Floating a) => Mat a
        mat = toMat [[1,1,0],[1,0,1],[0,1,1]]
    let qr :: (Integral a, RealFrac a, Floating a) => (Mat a, Mat a)
        qr = qrDecomp mat
    let q :: (Integral a, RealFrac a, Floating a) => Mat a
        q = fst qr
    let r :: (Integral a, RealFrac a, Floating a) => Mat a
        r = snd qr
    putStrLn ""
    putStrLn $ "Initial Matrix A ="
    print mat
    putStrLn "Matrix Q ="
    print q 
    putStrLn "Matrix R ="
    print r
    let multiplied :: (Integral a, RealFrac a, Floating a) => Mat a
        multiplied = mmMult r q
    if multiplied == mat
    then putStrLn "Q*R = A"
    else do 
        putStrLn "Q*R /= A"
        putStrLn $ "Q*R = " ++ (show multiplied)
        