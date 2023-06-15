-- |
-- Module      :  Bio.Reactamole.Real
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Functions for manipulating Reals in Reactamole.

module Bio.Reactamole.Real
  ( -- * Signal Functions
    -- ** Basic SFs
    negateSF
  , constRl
  , constAdd
  , constMult
  , sqrSF

    -- ** Comparators
  , isPosSF
  , ltSF
  , gtSF

    -- ** Real Function SFs
  , addSF
  , subSF
  , multSF
  , integrateSF
  , derivativeSF
  , tanhSF
  , signSF
  , recipSF
  , sqrtSF
  , logSF
  , sinSF
  , cosSF
  , expSF
  , coshSF
  , sinhSF
  , absSF
  , atanSF
  , rndSF
  , maxSF
  , minSF
  , lxhSF
  , selectSF
  , timeDilate
  , constTime
  , real2Bool
  ) where

import Bio.Reactamole.Arr
import Bio.Reactamole.Core
import Bio.Reactamole.Bool

--------------------------------------------------------------------------------

-- | Numerical negation for SFs.
negateSF :: SF Double Double
negateSF = arrTag (\(RealT eq) -> RealT (negateEq eq))

-- | Numerical addition for SFs.
addSF :: SF (Double, Double) Double
addSF = arrTag (\(PairT (RealT eq1) (RealT eq2)) -> RealT (addEq eq1 eq2))

-- | Numerical multiplication for SFs.
multSF :: SF (Double, Double) Double
multSF = arrTag (\(PairT (RealT eq1) (RealT eq2)) -> RealT (multEq eq1 eq2))

-- | Create an SF that ignores its input, producing the double d regardless.
constRl :: Double -> SF a Double
constRl d = constSF (Sg (RealT [Term d []]) [])

-- | Adds a constant to the signal
constAdd :: Double -> SF Double Double
constAdd d = constRl d &&& idSF >>> addSF

-- | Multiplies the signal by a constant
constMult :: Double -> SF Double Double
constMult d = constRl d &&& idSF >>> multSF

-- | Squares the signal
sqrSF :: SF Double Double
sqrSF = dupSF >>> multSF

-- | Numerical subtraction for SFs.
subSF :: SF (Double, Double) Double
subSF = second negateSF >>> addSF

--------------------------------------------------------------------------------

-- | Numerical integration for SFs from zero to infinity.
integrateSF :: Double -> SF Double Double
-- integrateSF c = SF $ \s@(Sg (RealT eq) _ _) -> addDualRail s eq c
integrateSF c = SF $ \(Sg (RealT eq) sys) ->
  let x = length sys
  in Sg (RealT [Term 1 [x]]) (sys ++ [ODE eq c Uncompiled (Just $ origin "integrateSF")])

--------------------------------------------------------------------------------

-- | Check if a double is positive.
isPosSF :: SF Double Bool
isPosSF = real2Bool

-- | Check if @f(t) > g(t)@ where @f@ and @g@ are the two input signals.
--
-- __WARNING__: Latency depends on how close @f(t)@ and @g(t)@ are to one
--              another and therefore may fail when @f(t) - g(t)@ is small.
gtSF :: SF (Double, Double) Bool
gtSF = subSF >>> isPosSF

-- | Check if @f(t) < g(t)@ where @f@ and @g@ are the two input signals.
--
-- __WARNING__: Latency depends on how close @f(t)@ and @g(t)@ are to one
--              another and therefore may fail when @f(t) - g(t)@ is small.
ltSF :: SF (Double, Double) Bool
ltSF = gtSF >>> notSF

--------------------------------------------------------------------------------

-- | Numerical differentiation for SFs
derivativeSF :: SF Double Double
derivativeSF = SF $ \(Sg (RealT x) sys) ->
  Sg (RealT (derivativeOfEq sys x)) sys

-- | Produces a SF that computes @f(x(t))@ where @x(t)@ is the input signal.
--
-- @f@ must be provided both as a real-valued function (to compute the initial
-- value from @x(0)@) and as a signal function computing its derivative
extendWith :: (Double -> Double) -> SF Double Double -> SF Double Double
extendWith f dfdt = SF $ \s@(Sg (RealT x) sys) ->
  runSF (dfdt >>> integrateSF (f (evalEq sys x))) s

-- | Produces a SF that computes @f(x(t))@ where @x(t)@ is the input signal.
--
-- @f@ must be provided both as a real-valued function (to compute the initial
-- value from @x(0)@) and as a signal function computing its derivative
--
-- __NOTE:__ 
extendWithLoop :: (Double -> Double) -> SF (Double, Double) Double -> SF Double Double
extendWithLoop f c = SF $ \s@(Sg (RealT x) sys) ->
  runSF (loop (c >>> integrateSF (f (evalEq sys x)) >>> dupSF)) s

-- | Numerical hyperbolic tangent function for SFs
tanhSF :: SF Double Double
tanhSF = extendWithLoop tanh dydt
  where dydt = derivativeSF *** (constRl 1 &&& sqrSF >>> subSF) >>> multSF

-- | Computes the reciprocal function
--
-- __WARNING__: Will fail if the input @f(t)@ is ever equal to zero
recipSF :: SF Double Double
recipSF = extendWithLoop (1/) dydt
  where dydt = derivativeSF *** (sqrSF >>> negateSF) >>> multSF

-- TODO: This is currently broken, but I can't figure out why
sqrtSF :: SF Double Double
sqrtSF = extendWithLoop sqrt dydt
  where dydt = derivativeSF *** (constMult 2 >>> recipSF) >>> multSF

-- | Computes the log function
--
-- __WARNING__: Will fail if the input @f(t)@ is ever less than or equal to zero
logSF :: SF Double Double
logSF = extendWith log dydt
  where dydt = derivativeSF &&& recipSF >>> multSF

sinCosSF :: SF Double (Double, Double)
sinCosSF = SF $ \s@(Sg (RealT eq) sys) ->
  let x0   = evalEq sys eq
      sinP = derivativeSF *** idSF     >>> multSF >>> integrateSF (sin x0)
      cosP = derivativeSF *** negateSF >>> multSF >>> integrateSF (cos x0)
  in runSF (loop (arrTag (\(PairT x (PairT sn cs)) -> PairT (PairT x cs) (PairT x sn))
             >>> (sinP *** cosP) >>> dupSF)) s

sinSF :: SF Double Double
sinSF = sinCosSF >>> proj1SF

cosSF :: SF Double Double
cosSF = sinCosSF >>> proj2SF

expSF :: SF Double Double
expSF = extendWithLoop exp dydt
  where dydt = first derivativeSF >>> multSF

coshSF :: SF Double Double
coshSF = expSF &&& (negateSF >>> expSF) >>> addSF >>> constMult (1/2)

sinhSF :: SF Double Double
sinhSF = expSF &&& (negateSF >>> expSF) >>> subSF >>> constMult (1/2)

absSF :: Double -> SF Double Double
absSF δ = constMult (1/δ) >>> tanhSF &&& idSF >>> multSF >>> constAdd δ

atanSF :: SF Double Double
atanSF = extendWith atan z
  where z    = extendWithLoop (\x -> 1/((x*x) + 1)) dzdt
        dzdt = (derivativeSF &&& idSF >>> multSF) *** sqrSF >>> multSF >>> constMult (-2)

--------------------------------------------------------------------------------

signSF :: Double -> Double -> SF Double Double
signSF μ λ = constMult (μ*λ) >>> tanhSF

ip1SF :: Double -> Double -> SF Double Double
ip1SF μ λ = constAdd (-1) >>> signSF μ λ >>> constAdd 1 >>> constMult (1/2)

rndSF :: SF Double Double
rndSF = idSF &&& (constMult c >>> sinSF >>> constMult (1/c)) >>> subSF
  where c = 2 * pi

maxSF :: Double -> SF (Double, Double) Double
maxSF δ = arrTag (\(PairT x y) -> PairT y x)
           >>> addSF &&& (subSF >>> absSF (2*δ))
           >>> addSF >>> constMult (1/2)

minSF :: Double -> SF (Double, Double) Double
minSF δ = addSF &&& maxSF δ >>> subSF

lxhSF :: Double -> SF (Double, Double) Double
lxhSF μ = (constAdd 1 >>> ip1SF μ 2) *** idSF >>> multSF

selectSF :: Double -> SF (Double, (Double, Double)) Double
selectSF μ = (second (subSF >>> negateSF) >>> lxhSF μ) &&& (proj2SF >>> proj1SF) >>> addSF

--------------------------------------------------------------------------------

-- | Time dilates the second SF by multiplying the output of the first SF
--   with the ODE of every variable used in the second SF
timeDilate :: SF a Double -> SF a b -> SF a b
timeDilate f g = f &&& g >>> SF (\(Sg (PairT (RealT x) y) sys) ->
    let yVars = depsTag y sys
        multIf i = let ODE eq ic cf src = sys !! i
                       eq' = if i `elem` yVars then multEq x eq else eq
                   in ODE eq' ic cf src
        sys' = map multIf [0..length sys - 1]
     in Sg y sys')

--------------------------------------------------------------------------------

real2Bool :: SF Double Bool
real2Bool = compileRl >>> SF (\(Sg (RealT [Term 1 [xp], Term (-1) [xn]]) sys) ->
  let k  = 30
      y  = length sys
      y' = y + 1
      rns = [ Rn [y,  y,  y'] [y,  y,  y ] (3*k)
            , Rn [y', y', y ] [y', y', y'] (3*k)
            , Rn [xn, y ] [xn, y'] k
            , Rn [xp, y'] [xp, y ] k ]
  in Sg (BoolT y y') (applyRns rns (sys ++ [ODE [] 0 Finalized src, ODE [] 1 Finalized src])))
  where
    src = Just $ origin "real2Bool"

--------------------------------------------------------------------------------

constTime :: SF a Double
constTime = constRl 1 >>> integrateSF 0
