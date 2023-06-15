-- |
-- Module      :  Bio.Reactamole.Export
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- A collection of useful functions for cleaning up, printing, and converting
-- between ODEs and SFs.

module Bio.Reactamole.Export
  (
    IVP (..)
  , CRN (..)
  , toIVP
  , toCRN
  , dualRail
  ) where

import Bio.Reactamole.Core
import Bio.Reactamole.Arr
import Data.List (delete, groupBy, sort, nub, intercalate)
import Data.Maybe (fromMaybe)
import qualified Data.Set as S

-- | An (IVP inSp outSp sys ic) object represents a complete initial value
--   problem where inSp helps interpret the input structure, outSp helps
--   interpret the output structure, sys represents the list of ODEs of the IVP,
--   and ic is the list of initial conditions
data IVP a b  = IVP (Tag a) (Tag b) [ODE]

data TgRxn = TgRxn { trSrc :: Maybe Origin, trRxn :: Reaction }
  deriving Eq

instance Show TgRxn where
  show (TgRxn Nothing rxn)    = show rxn
  show (TgRxn (Just src) rxn) = show rxn ++ " [" ++ intercalate "," (S.toList src) ++ "]"

instance Ord TgRxn where
  compare (TgRxn _ r1) (TgRxn _ r2) = compare r1 r2

-- | A (CRN inSp outSp rns ic) is similar except that it has a list of
--   reactions instead of a list of ODEs
data CRN a b = CRN (Tag a) (Tag b) [TgRxn] [Double]

-- | Creates a human-readable string describing the input tag
showIn :: Tag a -> String
showIn tg  = "INPUT:\n  "  ++ show tg ++ "\n"

-- | Creates a human-readable string describing the output tag
showOut :: Tag a -> String
showOut tg = "OUTPUT:\n  " ++ show tg ++ "\n"

-- | Creates a human-readable string describing the system of ODEs
--
--   NOTE: The input variables are hidden from the output since their ODEs
--         are undefined.
showSys :: [Variable] -> [ODE] -> String
showSys inputs sys = "ODEs:\n" ++ concat (zipWith displayODE [0..length sys-1] sys)
  where displayODE x (ODE eq ic cf s)
          | x `elem` inputs = ""
          | cf == Compiled = "  x" ++ show x ++ " = " ++ show ic ++ ", dx" ++ show x ++ "/dt = not available; use dualRail to finalize\n"
          | otherwise = "  x" ++ show x ++ " = " ++ show ic ++ ", dx" ++ show x ++ "/dt = " ++ show eq ++ src ++ "\n"
          where
            src = fromMaybe "" (fmap (\ns -> " [" ++ intercalate "," (S.toList ns) ++ "]") s)

-- | Creates a human-readable string describing the initial conditions of the systems
--
--   NOTE: The input variables are hidden from the output since their initial conditions
--         are undefined.
showICs :: [Variable] -> [Double] -> String
showICs inputs ic = "INITIAL CONDITIONS:\n" ++ concat (zipWith showIC [0..length ic - 1] ic)
  where showIC x x0 =
          if x `elem` inputs
            then ""
            else "  x" ++ show x ++ " = " ++ show x0 ++ "\n"

-- | Creates a human-readable string listing all the reactions given
showRns :: [TgRxn] -> String
showRns rns = "REACTIONS:\n" ++ intercalate "\n" (map (("  " ++) . show) rns) ++ "\n"

instance Show (IVP a b) where
  show (IVP ta tb sys) =
    intercalate "\n" [showIn ta, showOut tb, showSys inputs sys]
    where inputs = varsTag ta

instance Show (CRN a b) where
  show (CRN ta tb rns ic) =
    intercalate "\n" [showIn ta, showOut tb, showRns rns, showICs (varsTag ta) ic]

-- | Extract the initial value problem from the signal function.
toIVP :: HasDefault a => SF a b -> IVP a b
toIVP f = IVP inSp outSp (map normalizeODE sys)
  where s@(Sg inSp _) = runSF instSF (Sg getDefault [])
        f' = f >>> reduceSF (varsTag inSp)
        Sg outSp sys = runSF f' s

-- | Given a term in the ODE of a species x, returns the associated reaction
termToRn :: Variable -> Term -> Reaction
termToRn x (Term c ys) =
  if c > 0 then Rn ys (x:ys) c
           else Rn ys (delete x ys) (-c)

-- | Lifts @termToTgRxns@ to lists of terms
termsToRns :: Variable -> [Term] -> [Reaction]
termsToRns = map . termToRn

-- | Given a system of ODEs, returns the associated list of reactions
eqsToRns :: [(Maybe Origin, Equation)] -> [TgRxn]
eqsToRns sys = sort (map (foldr1 combineRns) (groupRns rns))
  where
    rns = concat $ zipWith
      (\v (src, terms) -> map (\rxn -> TgRxn src rxn) $ termsToRns v terms) [0..length sys - 1] sys

-- | Groups the given reactions by their *rates*, i.e., the left-hand
--   side of the reaction and the rate constant
groupRns :: [TgRxn] -> [[TgRxn]]
groupRns = groupBy f . sort
  where f (TgRxn _ (Rn r1 _ k1)) (TgRxn _ (Rn r2 _ k2)) = sort r1 == sort r2 && k1 == k2

-- | Combines two reactions with the same rate into a single reaction
--   by updating the products of the reaction
combineRns :: TgRxn -> TgRxn -> TgRxn
combineRns (TgRxn src1 (Rn r p1 k)) (TgRxn src2 (Rn _ p2 _)) = TgRxn (combineTags src1 src2) (Rn r p k)
  where tr  = varTally r
        tr' = map (\(x,n) -> (x,-n)) tr
        tp1 = varTally p1
        tp2 = varTally p2
        ne1 = addTallies tp1 tr'
        ne2 = addTallies tp2 tr'
        tp  = addTallies tr (addTallies ne1 ne2)
        p   = concatMap (\(x,n) -> replicate n x) tp
        combineTags t1 t2 = do
          s1 <- t1
          s2 <- t2
          pure $ S.union s1 s2

-- | Given two association lists mapping variables to integers,
--   produces a new association list that is the *sum* of the
--   two given lists.
addTallies :: [(Variable,Int)] -> [(Variable,Int)] -> [(Variable,Int)]
addTallies xs ys = map total (nub (map fst zs))
  where zs = xs ++ ys
        total z = (z,sum (map snd (filter ((==z) . fst) zs)))

finalizeODEs :: [Variable] -> [ODE] -> [ODE]
finalizeODEs inputs sys = map kernel [0..length sys - 1]
  where
    kernel v =
      let vp = if even v then v else v-1
          vn = vp + 1
          ODE eq _ cf src = sys !! vp
          extraTerm = if vp `elem` inputs then [] else [Term (-100) [vp, vn]]
          dvpdt = filter (\(Term c _) -> c > 0) eq ++ extraTerm
          dvndt = negateEq (filter (\(Term c _) -> c < 0) eq) ++ extraTerm
      in case (even v, cf) of
          (_    , Finalized  ) -> sys !! v
          (True , Compiled   ) -> ODE dvpdt (initial (sys !! v)) Finalized src
          (False, Compiled   ) -> ODE dvndt (initial (sys !! v)) Finalized src
          (_    , Uncompiled ) -> error "*** All ODEs must be compiled before this function is called"

dualRail :: IVP a b -> IVP a b
dualRail (IVP sa sb sys) = IVP sa' sb' (finalizeODEs inputs sys')
  where 
    sa' = expandTag sa
    sb' = expandTag sb
    inputs = varsTag sa'
    sys' = concat $ zipWith expandVar (map source sys) [0..length sys - 1]

    expandTag :: Tag a -> Tag a
    expandTag (RealT eq) = RealT (expandEq eq)
    expandTag (BoolT x x') = BoolT (shiftVar x) (shiftVar x')
    expandTag (PairT x y) = PairT (expandTag x) (expandTag y)
    expandTag (Tup3T x y z) = Tup3T (expandTag x) (expandTag y) (expandTag z)
    expandTag (Tup4T x y z w) = Tup4T (expandTag x) (expandTag y) (expandTag z) (expandTag w)
    expandTag (Tup5T x y z w u) = Tup5T (expandTag x) (expandTag y) (expandTag z) (expandTag w) (expandTag u)
    expandTag (CondT x y z) = CondT x (expandTag y) (expandTag z)
    expandTag tg = tg

    isCompiled :: Variable -> Bool
    isCompiled v = f == Compiled || f == Finalized where f = flag (sys !! v)

    shiftVar :: Variable -> Variable
    shiftVar v = v + sum (map (\i -> if isCompiled i then 0 else 1) [0..v-1])

    findDual :: Variable -> (Variable, Variable)
    findDual v = (v', v'+1) where v' = shiftVar v

    expandVar :: Maybe Origin -> Variable -> [ODE]
    expandVar src v
      | isCompiled v = let ODE dvdt v0 ic s = sys !! v in [ODE (expandEq dvdt) v0 ic s]
      | otherwise =
          let eq = expandEq (equation (sys !! v))
              vp0 = max (initial (sys !! v)) 0
              vn0 = max (-(initial (sys !! v))) 0
          in [ODE eq vp0 Compiled src, ODE [] vn0 Compiled src]
    expandEq :: Equation -> Equation
    expandEq = concatMap expandTerm

    expandTerm :: Term -> Equation
    expandTerm (Term c []) = [Term c []]
    expandTerm (Term c (x:xs))
      | isCompiled x = multEq [Term 1 [shiftVar x]] (expandTerm (Term c xs))
      | otherwise =
        let (xp, xn) = findDual x
        in multEq [Term 1 [xp], Term (-1) [xn]] (expandTerm (Term c xs))

-- | Converts the given signal function into an explicit CRN
toCRN :: HasDefault a => SF a b -> CRN a b
toCRN f = CRN ta tb (eqsToRns (map (\ode -> (source ode, equation ode)) sys)) (map initial sys)
  where IVP ta tb sys = dualRail (toIVP (f >>> instSF))

