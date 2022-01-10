-- |
-- Module      :  Bio.Reactamole.Examples
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Examples of uses for Reactamole.
--
--  * Including all examples presented in:
--    /Reactamole: Functional Reactive Molecular Programming/, by Kinge,
--    Lathrop, Osera and Rogers in /DNA Computing and Molecular Programming/ 27,
--    September 2021.

module Bio.Reactamole.Examples
  ( -- * Example CRNs
    sinCRN
  , expCRN
  , lowPass
  , bandPass
  , modulate
  , demodulate
  , srLatch
  , clock
  , unanimousCRN

    -- * Helpers
  , constMult
  , carrier
  , entangle
  , rectify
  , unanimous
  , inputSig
  , modulatedSig
  , demodulatedSig
  , withLowPassSig
  ) where

import Bio.Reactamole
import Bio.Reactamole.ArrChoice

-- | CRN for the sine function.
sinCRN :: CRN a Double
sinCRN = loop $ proj2 >>> negateCRN >>> integrate 1 >>> integrate 0 >>> dupCRN

-- | CRN for the exponential function.
expCRN :: CRN a Double
expCRN = loop (proj2 >>> integrate 1 >>> dupCRN)

-- | Multiply a 'Signal' by the given constant.
constMult :: Double -> CRN Double Double
constMult d = constRl d &&& idCRN >>> multCRN

-- | A lowpass filter.
lowPass :: Double -> Double -> CRN Double Double
lowPass a b = loop (constMult a *** constMult (-b) >>> addCRN >>> integrate 0 >>> dupCRN)

-- | A bandpass filter.
bandPass :: Double -> Double -> Double -> CRN Double Double
bandPass a b c = loop (first (constMult a)
  >>> second (constMult (-c) &&& (integrate 0 >>> constMult (-b)) >>> addCRN)
  >>> addCRN >>> integrate 0 >>> dupCRN)

-- | Generate a sine wave with the given frequency.
carrier :: Double -> CRN a Double
carrier w = loop (proj2 >>> constMult (-w) >>> integrate 1
        >>> constMult w >>> integrate 0 >>> dupCRN)

-- | Produce a pair @m(t)@ that satisfies
-- \( \frac{dm}{dt} = u(t)\cdot\sin(ft) - m(t) \).
modulate :: Double -> CRN Double Double
modulate w = loop (first (idCRN &&& carrier w >>> multCRN)
  >>> second negateCRN >>> addCRN >>> integrate 0 >>> dupCRN)

-- | Demodulate a signal that has been modulated on a carrier signal at
-- frequency @w@.
demodulate :: Double -> Double -> CRN Double Double
demodulate w q = bandPass (w/q) (w/q) (w*w) >>> rectify >>> lowPass (w/10) (w/10)

-- | Rectifies a signal so that only the positive part goes through
rectify :: CRN Double Double
rectify = posCRN &&& dupCRN >>> entangle >>> (idCRN ||| constRl 0)

--------------------------------------------------------------------------------

-- | An SR latch CRN.
srLatch :: CRN (Bool, Bool) (Bool, Bool)
srLatch = loop $
  arrSp (\(PairS (PairS s' r') (PairS q q'))
         -> PairS (PairS s' q') (PairS r' q))
  >>> (nandCRN *** nandCRN)
  >>> dupCRN

-- | A clock CRN.
clock :: CRN a Bool
clock = sinCRN >>> posCRN

-- | A Haskell function to determine if three inputs are the same.
--
-- For demonstration of lifting Boolean functions.
unanimous :: Bool -> Bool -> Bool -> Bool
unanimous x y z = if x then y && z else not (y || z)

-- | An CRN to determine if three inputs are the same.
--
-- Demonstrates the use of arr3Bl.
unanimousCRN :: CRN (Bool, Bool, Bool) Bool
unanimousCRN = arr3Bl unanimous

--------------------------------------------------------------------------------

exp1, exp2, exp3 :: CRN () Double
exp1 = constRl 0.5  &&& carrier 0.001  >>> multCRN
exp2 = constRl 0.25 &&& carrier 0.002  >>> multCRN
exp3 = constRl 0.3  &&& carrier 0.0014 >>> multCRN

inputSig :: CRN () Double
inputSig = constRl 1.1 &&& exp1 >>> addCRN
                       &&& exp2 >>> addCRN
                       &&& exp3 >>> addCRN
                &&& constRl 0.5 >>> multCRN

modulatedSig :: CRN () Double
modulatedSig = inputSig >>> modulate 0.75

demodulatedSig :: CRN () Double
demodulatedSig = modulatedSig >>> demodulate 0.75 0.25

withLowPassSig :: CRN () Double
withLowPassSig = demodulatedSig >>> lowPass 0.13 0.025
