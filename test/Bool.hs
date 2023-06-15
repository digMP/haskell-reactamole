module Bool where

import Test.Tasty
import Test.Tasty.HUnit

import Bio.Reactamole

truishConc :: [Double] -> Bool
truishConc [_, x] = x > eps
  where
    eps = 0.001
truishConc _      = error "(truishConc) given more than two concentrations"

signalBoolTestCase :: String -> Signal Bool -> Bool -> TestTree
signalBoolTestCase desc sg expected = testCase desc test
  where
    time = 2000.0
    test = do
      table <- calculateSignalTable time sg
      let results = last table
      truishConc results == expected @? "Expected: " ++ show expected ++ ", received: " ++ show results

allTestsUn :: String -> SF Bool Bool -> [Bool] -> TestTree
allTestsUn desc sf exps =
  testGroup desc $ zipWith3 mkCase allInputs [1..2] exps
    where
      allInputs = [trueSg, falseSg]
      mkCase initial n expected =
        signalBoolTestCase (show (n :: Integer))
                           (runSF (sf >>> instBl) initial)
                           expected

allTestsUnI :: String -> SF Bool Bool -> [Bool] -> TestTree
allTestsUnI desc sf exps =
  testGroup desc $ zipWith3 mkCase allInputs [1..2] exps
    where
      allInputs = [trueSg, falseSg]
      mkCase b n expected =
        signalBoolTestCase (show (n :: Integer))
                           (runSF (constSF b >>> instBl >>> sf) nullSg)
                           expected

allTestsBin :: String -> SF (Bool, Bool) Bool -> [Bool] -> TestTree
allTestsBin desc sf exps =
  testGroup desc $ zipWith3 mkCase allInputs [1..4] exps
  where
    allInputs =
      map (uncurry pairSg) [ (trueSg, trueSg)
                           , (trueSg, falseSg)
                           , (falseSg, trueSg)
                           , (falseSg, falseSg)
                           ]
    mkCase initial n expected =
      signalBoolTestCase (show (n :: Integer))
                         (runSF (sf >>> instBl) initial)
                         expected

allTestsBinI :: String -> SF (Bool, Bool) Bool -> [Bool] -> TestTree
allTestsBinI desc sf exps =
  testGroup desc $ zipWith3 mkCase allInputs [1..4] exps
  where
    allInputs = [ (trueSg, trueSg)
                , (trueSg, falseSg)
                , (falseSg, trueSg)
                , (falseSg, falseSg)
                ]
    mkCase (b1, b2) n expected =
      signalBoolTestCase (show (n :: Integer))
                         (runSF ((constSF b1 >>> instBl) &&& (constSF b2 >>> instBl) >>> sf) nullSg)
                         expected

--------------------------------------------------------------------------------

boolTests :: TestTree
boolTests = testGroup "Boolean Operators"
  [ unaryEarly
  , binaryEarly
  , unaryLate
  , binaryLate
  ]

unaryEarly :: TestTree
unaryEarly = testGroup "Unary (Early Instantiation)"
  [ allTestsUnI "idSF"  idSF  [True, False]
  , allTestsUnI "notSF" notSF [False, True]
  ]

binaryEarly :: TestTree
binaryEarly = testGroup "Binary (Early Instantiation)"
  [ allTestsBinI "nandSF" nandSF [False, True, True, True]
  , allTestsBinI "andSF"  andSF  [True, False, False, False]
  , allTestsBinI "orSF"   orSF   [True, True, True, False]
  , allTestsBinI "norSF"  norSF  [False, False, False, True]
  , allTestsBinI "xorSF"  xorSF  [False, True, True, False]
  ]

unaryLate :: TestTree
unaryLate = testGroup "Unary (Late Instantiation)"
  [ allTestsUn "idSF"  idSF  [True, False]
  , allTestsUn "notSF" notSF [False, True]
  ]

binaryLate :: TestTree
binaryLate = testGroup "Binary (Late Instantiation)"
  [ allTestsBin "nandSF" nandSF [False, True, True, True]
  , allTestsBin "andSF"  andSF  [True, False, False, False]
  , allTestsBin "orSF"   orSF   [True, True, True, False]
  , allTestsBin "norSF"  norSF  [False, False, False, True]
  , allTestsBin "xorSF"  xorSF  [False, True, True, False]
  ]
