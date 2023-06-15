-- |
-- Module      :  Bio.Reactamole.Arr
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Arrow definitions for SFs in Reactamole.
--
-- Sources:
--
--  * /Generalising Monads to Arrows/, by John Hughes,
--    /Science of Computer Programming/ 37, pp67-111, May 2000.
--
--  * /A New Notation for Arrows/, by Ross Paterson, in /ICFP 2001/,
--    Firenze, Italy, pp229-240.
--
--  * [Control.Arrow package.](https://hackage.haskell.org/package/base-4.15.0.0/docs/Control-Arrow.html)

module Bio.Reactamole.Arr
  ( -- * Arrows for SFs
    (<<<)
  , (>>>)
  , (***)
  , (&&&)
  , first
  , second
  , loop
  --, loop'

  , instSF

  --, derivativeOfTerm
  , derivativeOfEq
  , compileRl
  , defaultSg
  --, dropSg
  ) where

import Bio.Reactamole.Core

import qualified Data.Set as S

--------------------------------------------------------------------------------

-- | Right-to-left composition.
(<<<) :: SF b c -> SF a b -> SF a c
(SF f) <<< (SF g) = SF (f . g)

-- | Left-to-right composition.
(>>>) :: SF a b -> SF b c -> SF a c
f >>> g = g <<< f

-- | Split the input between the two argument arrows and combine their output.
--
-- Note that this is in general not a functor.
(***) :: SF a b -> SF a' b' -> SF (a, a') (b, b')
(SF f) *** (SF g) = SF $ \s@(Sg _ sys) ->
  let n = length sys
      (Sg spf sysf) = (dropSg n . f . runSF proj1SF) s
      k = length sysf
      (Sg spg sysg) = (shiftSg n k . dropSg n . g . runSF proj2SF) s
    in Sg (PairT spf spg) (sys ++ sysf ++ sysg)

-- | Fanout: send the input to both argument arrows and combine their output.
(&&&) :: SF a b -> SF a c -> SF a (b, c)
f &&& g = dupSF >>> (f *** g)

-- | Send the first component of the input through the argument arrow, and copy
-- the rest unchanged to the output.
first :: SF a b -> SF (a, c) (b, c)
first f = f *** idSF

-- | A mirror image of 'first'.
second :: SF a b -> SF (c, a) (c, b)
second f = idSF *** f

--------------------------------------------------------------------------------

-- | Loop the second parameter, c, back into itself.
--
-- The c in the resulting SF is unused and unmodified by the new SF
-- (essentially the identity).
loop' :: SF (a,c) (b,c) -> SF (a,c) (b,c)
loop' f = SF $ \s@(Sg (PairT _ y) _) ->
  let s'@(Sg (PairT _ y') _) = runSF (second instSF) s
      s''@(Sg (PairT _ y'') _) =  runSF (f >>> second instSF) s'
      Sg (PairT x _) sys = rewireVars s'' $ zip (varsTag y'') (varsTag y')
  in Sg (PairT x y) sys

-- | Feed the second component back in on itself.
loop :: HasDefault c => SF (a, c) (b, c) -> SF a b
loop f = SF $ \(Sg x sys) ->
  let SF g = loop' f
      Sg (PairT y _) sys' = g (Sg (PairT x getDefault) sys)
  in Sg y sys'

--------------------------------------------------------------------------------

-- | Instantiate all /uninstantiated/ variables.
instSF :: SF a a
instSF = SF $ \s@(Sg tg _) -> case tg of
  NullT     -> s
  TrueT     -> runSF instBl s
  FalseT    -> runSF instBl s
  BoolT _ _ -> runSF instBl s
  RealT _   -> runSF instRl s
  PairT _ _ -> runSF (instSF *** instSF) s
  Tup3T {}  -> runSF instTup3 s
  Tup4T {}  -> runSF instTup4 s
  Tup5T {}  -> runSF instTup5 s
  CondT {}  -> runSF instCond s

-- | Instantiate 'BoolT' variables.
instBl :: SF Bool Bool
instBl = SF $ \(Sg tg sys) ->
  let n = length sys
  in case tg of
    TrueT     -> Sg (BoolT n (n+1)) (sys ++ [ODE [] 1 Finalized src, ODE [] 0 Finalized src])
    FalseT    -> Sg (BoolT n (n+1)) (sys ++ [ODE [] 0 Finalized src, ODE [] 1 Finalized src])
    BoolT _ _ -> Sg tg sys
  where
    src = Just $ origin "instBl"

-- | Instantiate 'RealT' signal to be a concrete single-rail variable
instRl :: SF Double Double
instRl = SF $ \(Sg (RealT eq) sys) ->
  case eq of
    [Term 1 [_]] -> Sg (RealT eq) sys
    _ -> Sg (RealT [Term 1 [length sys]])
            (sys ++ [ODE (derivativeOfEq sys eq) (evalEq sys eq) Uncompiled src])
  where
    src = Just $ origin "instRl"

-- | Instaniates a 'RealT' signal to be a dual-rail variable pair
compileRl :: SF Double Double
compileRl = SF $ \(Sg (RealT eq) sys) ->
  let x    = length sys
      x'   = x + 1
      x0   = evalEq sys eq
      tg   = RealT [Term 1 [x], Term (-1) [x']]
      sys' = sys ++ [ ODE (derivativeOfEq sys eq) (max x0 0) Compiled src
                    , ODE [] (max (-x0) 0) Compiled src
                    ]
  in Sg tg sys'
  where
    src = Just $ origin "compileRl"


-- | Instantiate 'Tup3T' variables.
instTup3 :: SF (a,b,c) (a,b,c)
instTup3 = arrTag (\(Tup3T x y z) -> PairT x (PairT y z))
             >>> (instSF *** instSF)
             >>> arrTag (\(PairT x (PairT y z)) -> Tup3T x y z)

-- | Instantiate 'Tup4T' variables.
instTup4 :: SF (a,b,c,d) (a,b,c,d)
instTup4 = arrTag (\(Tup4T x y z w) -> PairT x (Tup3T y z w))
             >>> (instSF *** instSF)
             >>> arrTag (\(PairT x (Tup3T y z w)) -> Tup4T x y z w)

-- | Instantiate 'Tup5T' variables.
instTup5 :: SF (a,b,c,d,e) (a,b,c,d,e)
instTup5 = arrTag (\(Tup5T x y z w u) -> PairT x (Tup4T y z w u))
             >>> (instSF *** instSF)
             >>> arrTag (\(PairT x (Tup4T y z w u)) -> Tup5T x y z w u)

-- | Instantiate 'CondT' variables.
instCond :: SF (Either a b) (Either a b)
instCond = arrTag (\(CondT x y z) -> PairT x (PairT y z))
             >>> (instSF *** instSF)
             >>> arrTag (\(PairT x (PairT y z)) -> CondT x y z)

--------------------------------------------------------------------------------

-- | Drop the first n variables of a signal.
--
-- Remove the first n _equations_ and the first n _initial conditions_ from the
-- initial value problem
--
-- __WARNING__: This does not re-index the variables! Use with caution!
--              Consider using with 'shiftSg' to re-index them.
dropSg :: Int -> Signal a -> Signal a
dropSg n (Sg tg sys) = Sg tg (drop n sys)

--------------------------------------------------------------------------------

-- | Instantiates the output signal of the SF by giving it a default (empty)
--   input signal
defaultSg :: HasDefault a => SF a b -> Signal b
defaultSg f = runSF (f >>> instSF >>> reduceSF []) s
  where s = runSF instSF (Sg getDefault [])
