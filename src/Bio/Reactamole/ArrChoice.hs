-- |
-- Module      :  Bio.Reactamole.ArrChoice
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Arrow definitions for SF conditionals and 'Bio.Reactamole.Core.Either'
-- values.
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

module Bio.Reactamole.ArrChoice
  ( -- * Conditional SFs
    entangle
  , split
  , (+++)
  , (|||)
  , left
  , right
  ) where

import Bio.Reactamole.Bool
import Bio.Reactamole.Real
import Bio.Reactamole.Arr
import Bio.Reactamole.Core

--------------------------------------------------------------------------------

split :: SF a Bool -> SF a (Either a a)
split f = f &&& (idSF &&& idSF) >>> entangle

-- | The entangle signal function.
entangle :: SF (Bool, (a, b)) (Either a b)
entangle = arrTag (\(PairT x (PairT y z)) -> CondT x y z)

--------------------------------------------------------------------------------

-- | Split the input between the two argument 'SF's, retagging and merging
-- their outputs.
--
-- Note that this is in general not a functor.
(+++) :: SF a b -> SF c d -> SF (Either a c) (Either b d)
f +++ g = arrTag (\(CondT x y z) -> PairT x (PairT y z))
  >>> (idSF *** (f *** g))
  >>> arrTag (\(PairT x (PairT y z)) -> CondT x y z)

-- | Fanin: split the input between the two argument arrows and merge their
-- outputs.
(|||) :: SF a c -> SF b c -> SF (Either a b) c
f ||| g = (f +++ g) >>> mergeSF

-- | Feed marked inputs through the argument arrow, passing the rest through
-- unchanged to the output.
left :: SF a b -> SF (Either a c) (Either b c)
left f = f +++ idSF

-- | A mirror image of left.
right :: SF a b -> SF (Either c a) (Either c b)
right f = idSF +++ f

--------------------------------------------------------------------------------

-- | Merge a conditional bool signal.
--
-- Implmentation uses the equivalence of the following two statements:
--   * @w = if x then y else z@
--   * @w = (x && y) || (!x && z)@
mergeBool :: SF (Either Bool Bool) Bool
mergeBool = arrTag (\(CondT x y z) ->
  PairT (PairT x y) (PairT x z))
  >>> (andSF *** (first notSF >>> andSF))
  >>> orSF

-- | Merge a conditional real signal.
--
-- Uses the common GPAC switching function that is accurate to within a
-- concentration of 0.05, as long as the Boolean has a value between (0.5, 1)
-- or (-1, -0.5)
mergeReal :: SF (Either Double Double) Double
mergeReal = arrTag (\(CondT (BoolT x x') y z) ->
  PairT (RealT [Term 1 [x], Term (-1) [x']]) (PairT y z))
  >>> selectSF 3

-- | Merge a conditional pair signal.
mergePair :: SF (Either (a, b) (a, b)) (a, b)
mergePair = arrTag (\(CondT x (PairT ya yb) (PairT za zb)) ->
  PairT (CondT x ya za) (CondT x yb zb))
  >>> (mergeSF *** mergeSF)

-- | Merge a conditional tuple 3 signal.
mergeTup3 :: SF (Either (a, b, c) (a, b, c)) (a, b, c)
mergeTup3 = arrTag (\(CondT x (Tup3T y z w) (Tup3T y' z' w')) ->
  PairT (CondT x y y') (CondT x (PairT z w) (PairT z' w')))
  >>> (mergeSF *** mergePair)
  >>> arrTag (\(PairT x (PairT y z)) -> Tup3T x y z)

-- | Merge a conditional tuple 4 signal.
mergeTup4 :: SF (Either (a, b, c, d) (a, b, c, d)) (a, b, c, d)
mergeTup4 = arrTag (\(CondT x (Tup4T y z w u) (Tup4T y' z' w' u')) ->
  PairT (CondT x y y') (CondT x (Tup3T z w u) (Tup3T z' w' u')))
  >>> (mergeSF *** mergeTup3)
  >>> arrTag (\(PairT x (Tup3T y z w)) -> Tup4T x y z w)

-- | Merge a conditional tuple 5 signal.
mergeTup5 :: SF (Either (a, b, c, d, e) (a, b, c, d, e)) (a, b, c, d, e)
mergeTup5 = arrTag (\(CondT x (Tup5T y z w u v) (Tup5T y' z' w' u' v')) ->
  PairT (CondT x y y') (CondT x (Tup4T z w u v) (Tup4T z' w' u' v')))
  >>> (mergeSF *** mergeTup4)
  >>> arrTag (\(PairT x (Tup4T y z w u)) -> Tup5T x y z w u)

-- | Merge a nested conditional signal.
--
-- Cosnider the following conditional:
--  * @if x then (if y then y' else y'')
--          else (if z then z' else z'')@
--
-- The above conditional is converted to: (@if w then w' else w''@) where
--  * @w@   is the merged Bool:      (@if x then y   else z  @)
--  * @w'@  is the merged Tag a: (@if x then y'  else z' @)
--  * @w''@ is the merged Tag b: (@if x then y'' else z''@)
mergeCond :: SF (Either (Either a b) (Either a b)) (Either a b)
mergeCond = arrTag (\(CondT x (CondT y y' y'') (CondT z z' z'')) ->
  PairT (CondT x y z) (PairT (CondT x y' z')
                             (CondT x y'' z'')))
  >>> (mergeBool *** (mergeSF *** mergeSF))
  >>> arrTag (\(PairT x (PairT y z)) -> CondT x y z)

-- | Merge a conditional signal.
mergeSF :: SF (Either c c) c
mergeSF = SF $ \s@(Sg (CondT x y z) sys) -> case (x,y,z) of
  (TrueT , _        , _) -> Sg y sys
  (FalseT, _        , _) -> Sg z sys
  (_     , NullT    , _) -> Sg NullT []
  (_     , TrueT    , _) -> runSF mergeBool s
  (_     , FalseT   , _) -> runSF mergeBool s
  (_     , BoolT {} , _) -> runSF mergeBool s
  (_     , RealT _  , _) -> runSF mergeReal s
  (_     , PairT {} , _) -> runSF mergePair s
  (_     , Tup3T {} , _) -> runSF mergeTup3 s
  (_     , Tup4T {} , _) -> runSF mergeTup4 s
  (_     , Tup5T {} , _) -> runSF mergeTup5 s
  (_     , CondT {} , _) -> runSF mergeCond s
