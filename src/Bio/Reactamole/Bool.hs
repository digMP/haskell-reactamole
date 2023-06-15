-- |
-- Module      :  Bio.Reactamole.Bool
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Tools for Bools: Functions for lifting and manipulating Bools in
-- Reactamole.
--
-- Source: NAND gate construction from /Robust Chemical Circuits/, by Ellis,
-- Klinge and Lathrop, /Biosystems/, 2019.

module Bio.Reactamole.Bool
  ( -- * Signal Functions
    -- ** Constant Boolean
    constBl

    -- ** Logic Gates
  , notSF
  , nandSF
  , andSF
  , orSF
  , norSF
  , xorSF
  , xnorSF

    -- ** Lift SFs
  , arr1Bl
  , arr2Bl
  , arr3Bl
  , arr4Bl
  , arr5Bl
  --, arrNOT
  --, arrAND
  --, arrOR
  ) where

import Bio.Reactamole.Arr
import Bio.Reactamole.Core

import qualified Data.Set as S

-- | SF that produces True (/'instantiate'd/) constantly.
constBl :: Bool -> SF a Bool
constBl b = constSF (Sg b' [])
  where b' = if b then TrueT else FalseT

--------------------------------------------------------------------------------

-- | SF for logical NOT.
notSF :: SF Bool Bool
notSF = arrTag $ \tg -> case tg of
  TrueT     -> FalseT
  FalseT    -> TrueT
  BoolT x y -> BoolT y x

--------------------------------------------------------------------------------

-- | SF for logical NAND.
nandSF :: SF (Bool, Bool) Bool
nandSF = SF $ \(Sg tg sys) -> case tg of
  PairT TrueT          TrueT          -> Sg FalseT         sys
  PairT FalseT         _              -> Sg TrueT          sys
  PairT _              FalseT         -> Sg TrueT          sys
  PairT (BoolT x1 x1') TrueT          -> Sg (BoolT x1' x1) sys
  PairT TrueT          (BoolT x2 x2') -> Sg (BoolT x2' x2) sys
  PairT (BoolT x1 x1') (BoolT x2 x2') ->
    let k  = 30
        y  = length sys
        y' = y + 1
        rns = [ Rn [y,   y,  y'] [y,   y,  y ] (3*k)
              , Rn [y',  y', y ] [y',  y', y'] (3*k)
              , Rn [x1,  x2, y ] [x1,  x2, y'] k
              , Rn [x1', y'    ] [x1', y     ] k
              , Rn [x2', y'    ] [x2', y     ] k     ]
    in Sg (BoolT y y') (applyRns rns (sys ++ [ODE [] 0 Finalized src, ODE [] 1 Finalized src])) 
  where
    src = Just $ origin "nandSF"

--------------------------------------------------------------------------------

-- | SF for logical AND.
andSF :: SF (Bool, Bool) Bool
andSF = nandSF >>> notSF

-- | SF for logical OR.
orSF :: SF (Bool, Bool) Bool
orSF = (notSF *** notSF) >>> nandSF

-- | SF for logical NOR.
norSF :: SF (Bool, Bool) Bool
norSF = orSF >>> notSF

-- | SF for logical XOR.
--
-- Note that this uses three nandSFs under the hood. All other gates except
-- xnorSF use only one.
xorSF :: SF (Bool, Bool) Bool
xorSF = (nandSF &&& orSF) >>> andSF

-- | SF for logical XNOR.
--
-- Note that this uses three nandSFs under the hood. All other gates except
-- xorSF use only one.
xnorSF :: SF (Bool, Bool) Bool
xnorSF = xorSF >>> notSF

--------------------------------------------------------------------------------

-- | Produce either an identical or inverse signal depending on the Boolean.
--
-- Lift helper function.
arrNOT :: Bool -> SF Bool Bool
arrNOT x = if x then idSF else notSF

-- | Compile a series of AND and NOT gates into an SF, to be fed into OR gates
--  that will simulate a Boolean function.
--
-- Lift helper function.
arrAND :: [Bool] -> Tag a -> SF a Bool
arrAND [] _ = constBl False
arrAND [x] (BoolT _ _) = arrNOT x
arrAND [_] _ = error "Tag must be a Bool or Tuple of Bools, and no longer than the Bool list"
arrAND (x:xs) (PairT (BoolT _ _) b@(BoolT _ _)) =
  andSF <<< (arrNOT x *** arrAND xs b)
arrAND (x:xs) (Tup3T (BoolT _ _) b@(BoolT _ _) c@(BoolT _ _)) =
  andSF <<< (arrNOT x *** arrAND xs (PairT b c)) <<< tup3ToPairSF
arrAND (x:xs) (Tup4T (BoolT _ _) b@(BoolT _ _) c@(BoolT _ _) d@(BoolT _ _)) =
  andSF <<< (arrNOT x *** arrAND xs (Tup3T b c d)) <<< tup4ToPairSF
arrAND (x:xs) (Tup5T (BoolT _ _)
               b@(BoolT _ _) c@(BoolT _ _) d@(BoolT _ _) e@(BoolT _ _)) =
  andSF <<< (arrNOT x *** arrAND xs (Tup4T b c d e)) <<< tup5ToPairSF
arrAND (_:_) _ = error "Bool list cannot be longer than Tag; Tag must be Bool or a tuple of Bools"

-- | Compile a series of OR, AND, and NOT gates into an SF to simulate a Boolean
--  function.
--
-- Lift helper function.
arrOR :: [[Bool]] -> Tag a -> SF a Bool
arrOR [] s = arrAND [] s
arrOR [x] s = arrAND x s
arrOR (x:xs) s = orSF <<< (arrAND x s *** arrOR xs s) <<< dupSF

--------------------------------------------------------------------------------

-- | Lift a 1-in 1-out Boolean Haskell function into an 'SF'.
arr1Bl :: (Bool -> Bool) -> SF Bool Bool
arr1Bl f = SF $ \s@(Sg tg _) ->
     let bl = [True, False]
         matrix = [(f x, [x]) | x <- bl]
         mFilt  = map snd . filter ((==True) . fst)
         SF sg = arrOR (mFilt matrix) tg
     in sg s

-- | Lift a 2-in 1-out Boolean Haskell function into an 'SF'.
arr2Bl :: (Bool -> Bool -> Bool) -> SF (Bool, Bool) Bool
arr2Bl f = SF $ \s@(Sg tg _) ->
  let bl = [True, False]
      matrix = [(f x y, [x, y]) | x <- bl, y <- bl]
      mFilt  = map snd . filter ((==True) . fst)
      SF sg = arrOR (mFilt matrix) tg
  in sg s

-- | Lift a 3-in 1-out Boolean Haskell function into an 'SF'.
arr3Bl :: (Bool -> Bool -> Bool -> Bool) -> SF (Bool, Bool, Bool) Bool
arr3Bl f = SF $ \s@(Sg tg _) ->
  let bl = [True, False]
      matrix = [(f x y z, [x, y, z]) | x <- bl, y <- bl, z <- bl]
      mFilt  = map snd . filter ((==True) . fst)
      SF sg = arrOR (mFilt matrix) tg
  in sg s
      
-- | Lift a 4-in 1-out Boolean Haskell function into an 'SF'.
arr4Bl :: (Bool -> Bool -> Bool -> Bool -> Bool) ->
          SF (Bool, Bool, Bool, Bool) Bool
arr4Bl f = SF $ \s@(Sg tg _) ->
  let bl = [True, False]
      matrix = [(f x y z w, [x, y, z, w]) |
                x <- bl, y <- bl, z <- bl, w <- bl]
      mFilt  = map snd . filter ((==True) . fst)
      SF sg = arrOR (mFilt matrix) tg
  in sg s
  
-- | Lift a 5-in 1-out Boolean Haskell function into an 'SF'.
arr5Bl :: (Bool -> Bool -> Bool -> Bool -> Bool -> Bool) ->
          SF (Bool, Bool, Bool, Bool, Bool) Bool
arr5Bl f = SF $ \s@(Sg tg _) ->
  let bl = [True, False]
      matrix = [(f x y z w v, [x, y, z, w, v]) |
                x <- bl, y <- bl, z <- bl, w <- bl, v <- bl]
      mFilt  = map snd . filter ((==True) . fst)
      SF sg = arrOR (mFilt matrix) tg
  in sg s
