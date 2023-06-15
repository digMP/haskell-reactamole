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
  ( -- * Example SFs
    lowPass
  , bandPass
  , modulate
  , demodulate
  , demodulate'
  , srLatch
  , clock
  , unanimousSF

    -- * Helpers
  , carrier
  , carrier'
  , entangle
  , rectify
  , rectify'
  , unanimous
  , inputSig
  , modulatedSig
  , demodulatedSig
  , withLowPassSig
  ) where

import Bio.Reactamole
import Bio.Reactamole.ArrChoice

-- | Rectifies a signal so that only the positive part goes through
rectify :: SF Double Double
rectify = isPosSF &&& dupSF >>> entangle >>> (idSF ||| constRl 0)

rectify' :: SF Double Double
rectify' = constRl 0 &&& idSF >>> maxSF 0.1

-- | Generate a sine wave with the given frequency.
carrier :: Double -> SF a Double
carrier w = loop (proj2SF >>> constMult (-w) >>> integrateSF 1
                          >>> constMult w    >>> integrateSF 0
                          >>> dupSF)

carrier' :: Double -> SF a Double
carrier' w = timeDilate (constRl w) (constTime >>> sinSF)

-- | A lowpass filter.
lowPass :: Double -> Double -> SF Double Double
lowPass a b = loop (constMult a *** constMult (-b) >>> addSF >>>
                    integrateSF 0 >>> dupSF)

-- | A bandpass filter.
bandPass :: Double -> Double -> Double -> SF Double Double
bandPass a b c = loop (first (constMult a)
  >>> second (constMult (-c) &&& (integrateSF 0 >>> constMult (-b)) >>>
              addSF)
  >>> addSF >>> integrateSF 0 >>> dupSF)

-- | Produce a pair @m(t)@ that satisfies
-- \( \frac{dm}{dt} = u(t)\cdot\sin(ft) - m(t) \).
modulate :: Double -> SF Double Double
modulate w = loop (first (idSF &&& carrier w >>> multSF)
  >>> second negateSF >>> addSF >>> integrateSF 0 >>> dupSF)

-- | Demodulate a signal that has been modulated on a carrier signal at
-- frequency @w@.
demodulate :: Double -> Double -> SF Double Double
demodulate w q = bandPass (w/q) (w/q) (w*w) >>> rectify >>> lowPass (w/10) (w/10)

demodulate' :: Double -> Double -> SF Double Double
demodulate' w q = bandPass (w/q) (w/q) (w*w) >>> rectify' >>> lowPass (w/10) (w/10)

--------------------------------------------------------------------------------

-- | An SR latch SF.
srLatch :: SF (Bool, Bool) (Bool, Bool)
srLatch = loop $
  arrTag (\(PairT (PairT s' r') (PairT q q'))
         -> PairT (PairT s' q') (PairT r' q))
  >>> (nandSF *** nandSF)
  >>> dupSF

-- | A clock SF.
clock :: SF a Bool
clock = constTime >>> sinSF >>> real2Bool

-- | A Haskell function to determine if three inputs are the same.
--
-- For demonstration of lifting Boolean functions.
unanimous :: Bool -> Bool -> Bool -> Bool
unanimous x y z = if x then y && z else not (y || z)

-- | An SF to determine if three inputs are the same.
--
-- Demonstrates the use of arr3Bl.
unanimousSF :: SF (Bool, Bool, Bool) Bool
unanimousSF = arr3Bl unanimous

--------------------------------------------------------------------------------

exp1, exp2, exp3 :: SF () Double
exp1 = constRl 0.5  &&& carrier 0.001  >>> multSF
exp2 = constRl 0.25 &&& carrier 0.002  >>> multSF
exp3 = constRl 0.3  &&& carrier 0.0014 >>> multSF

inputSig :: SF () Double
inputSig = constRl 1.1 &&& exp1 >>> addSF
                       &&& exp2 >>> addSF
                       &&& exp3 >>> addSF
                &&& constRl 0.5 >>> multSF

modulatedSig :: SF () Double
modulatedSig = inputSig >>> modulate 0.75 >>> fromSource "modulatedSig"

demodulatedSig :: SF () Double
demodulatedSig = modulatedSig >>> demodulate 0.75 0.25

withLowPassSig :: SF () Double
withLowPassSig = demodulatedSig >>> lowPass 0.13 0.025
