{-# LANGUAGE DeriveFunctor, ScopedTypeVariables #-}

-- |
-- Module      :  Bio.Reactamole.Core
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Core Reactamole language structure based on a representation of chemical
-- reaction networks (SFs) as ordinary differential equations (ODEs).

module Bio.Reactamole.Core
  (  -- * Signal Functions
    SF (..)
  , reduceSF
  , fromSource

    -- ** Basic SFs
  , idSF
  , proj1SF
  , proj2SF
  , constSF
  , dupSF

    -- ** Tuple/Pair conversions
  , tup3ToPairSF
  , tup4ToPairSF
  , tup5ToPairSF

    -- * Tag & Signals
    -- ** Construction
  , Tag (..)
  , Signal (..)
  , HasDefault (..)
  , arrTag

    -- ** Manipulation
  , varsTag
  , shiftSg
  
    -- * Ordinary Differential Equations
    -- ** Construction
  , Variable
  , Term (..)
  , Equation
  , CompileFlag (..)
  , Origin
  , origin
  , ODE (..)

    -- ** Manipulation
  , showEq
  , constEq
  , addEq
  , multEq
  , evalEq
  , negateEq
  , varsEq
  , derivativeOfEq
  , derivativeOfTerm
  , normalizeEq
  , normalizeODE
  , depsTag
  , rewireVars
  , varTally

    -- * Reactions
  , Reaction (..)
  , applyRns
  ) where

import Data.List (intercalate, nub, sort)
import qualified Data.Set as S

--------------------------------------------------------------------------------

-- | Variable of an ODE.
--
-- Acts as an index representing a unique variable name.
type Variable = Int

-- | Term of an ODE.
--
-- For example, @Term 5 [1, 1, 2, 1]@ represents \( \frac{dx}{dt} = 5x^3y \),
-- where \( x \) is represented by 'Variable' @1@ and \( y \) by 'Variable' @2@.
data Term = Term { coef :: Double -- ^Coefficient
                 , vars :: [Variable] -- ^Variables
                 }
            deriving (Eq)

instance Show Term where
  show (Term c vs) = sgn ++ term
    where sgn = if c < 0 then "-" else "+"
          term = intercalate "*x" (show (abs c) : map showVar (varTally vs))
          showVar (a, b) =
            if b == 1 then show a else show a ++ "^" ++ show b

instance Ord Term where
  compare (Term _ v1) (Term _ v2) = compare (sort v1) (sort v2)

-- | Tally a list of variables and remove duplicates.
varTally :: [Variable] -> [(Variable, Int)]
varTally vs = [(x, sum [1 | y <- vs, y == x]) | x <- nub vs]

--------------------------------------------------------------------------------

-- | Ordinary differential equation (ODE).
--
-- For example, @Equation [('Term' 3 [1, 2]), ('Term' -1 [1, 1, 3])]@ represents
-- \( dx\/dt = +3*x*y -x*x*z \), where \( x \) is represented by 'Vairable' @1@,
-- \( y \) by 'Variable' @2@, and \( z \) by 'Variable' @3@.
type Equation = [Term]

-- | Turn a polynomial equation into a human readable form.
showEq :: Equation -> String
showEq ts = unwords $ map show ts

-- | Returns a constant equation
constEq :: Double -> Equation
constEq d = [Term d []]

-- | Add two polynomial equations together.
addEq :: Equation -> Equation -> Equation
addEq = (++)

-- | Multiply two polynomial equations together.
multEq :: Equation -> Equation -> Equation
multEq eq1 eq2 = [multTerm t1 t2 | t1 <- eq1, t2 <- eq2]
  where multTerm (Term c1 v1) (Term c2 v2) = Term (c1*c2) (v1 ++ v2)

-- | Negate a polynomial equation.
negateEq :: Equation -> Equation
negateEq eq = multEq eq [Term (-1) []]

-- | Evaluate an equation with initial conditions.
evalEq  :: [ODE] -> Equation -> Double
evalEq sys eq = kernel eq 0
  where kernel [] total = total
        kernel ((Term c vs):ts) total =
          let val = c * product [initial (sys !! v) | v <- vs]
          in kernel ts (total + val)

-- | Return derivative of the 'Term'.
derivativeOfTerm :: [ODE] -> Term -> Equation
derivativeOfTerm _ (Term _ []) = []
derivativeOfTerm sys (Term c (x:xs)) =
  let x'  = equation(sys !! x)
      xs' = derivativeOfTerm sys (Term 1 xs)
  in addEq (multEq x' [Term c xs]) (multEq [Term c [x]] xs')

-- | Return derivative of the 'Equation'.
derivativeOfEq :: [ODE] -> Equation -> Equation
derivativeOfEq sys = concatMap (derivativeOfTerm sys)


--------------------------------------------------------------------------------

-- | Type tg definition.
data Tag a where
  NullT  :: Tag ()   -- ^ Empty type
  TrueT  :: Tag Bool -- ^ Boolean True type /(uninstantiated)/
  FalseT :: Tag Bool -- ^ Boolean False type /(uninstantiated)/
  BoolT  :: Variable  -> Variable  -> Tag Bool
  RealT  :: Equation  -> Tag Double
  PairT  :: Tag a -> Tag b -> Tag (a, b)
  Tup3T  :: Tag a -> Tag b -> Tag c -> Tag (a, b, c)
  Tup4T  :: Tag a -> Tag b -> Tag c -> Tag d -> Tag (a, b, c, d)
  Tup5T  :: Tag a -> Tag b -> Tag c -> Tag d -> Tag e -> Tag (a, b, c, d, e)
  CondT  :: Tag Bool -> Tag a -> Tag b -> Tag (Either a b) -- ^ Conditional type

--------------------------------------------------------------------------------

-- | Default value for a 'Tag'.
class HasDefault a where
  getDefault :: Tag a

instance HasDefault () where
  getDefault = NullT

instance HasDefault Bool where
  getDefault = FalseT

instance HasDefault Double where
  getDefault = RealT []

instance (HasDefault a, HasDefault b) => HasDefault (a,b) where
  getDefault = PairT (getDefault :: Tag a) (getDefault :: Tag b)

instance (HasDefault a, HasDefault b, HasDefault c)
  => HasDefault (a,b,c) where
  getDefault = Tup3T (getDefault :: Tag a)
                     (getDefault :: Tag b)
                     (getDefault :: Tag c)

instance (HasDefault a, HasDefault b, HasDefault c, HasDefault d)
  => HasDefault (a,b,c,d) where
  getDefault = Tup4T (getDefault :: Tag a)
                     (getDefault :: Tag b)
                     (getDefault :: Tag c)
                     (getDefault :: Tag d)

instance (HasDefault a, HasDefault b, HasDefault c, HasDefault d, HasDefault e)
  => HasDefault (a,b,c,d,e) where
  getDefault = Tup5T (getDefault :: Tag a)
                     (getDefault :: Tag b)
                     (getDefault :: Tag c)
                     (getDefault :: Tag d)
                     (getDefault :: Tag e)

instance (HasDefault a, HasDefault b) => HasDefault (Either a b) where
  getDefault = CondT (getDefault :: Tag Bool)
                     (getDefault :: Tag a)
                     (getDefault :: Tag b)

--------------------------------------------------------------------------------

instance Show (Tag a) where
  show tg =
    case tg of
      NullT           -> "Null"
      TrueT           -> "True"
      FalseT          -> "False"
      BoolT x y       -> concat ["Bool(x", show x, ", x", show y, ")"]
      RealT eq        -> concat ["Real(", showEq (normalizeEq eq), ")"]
      CondT x y z     -> concat ["Cond(", show x, show y, show z, ")"]
      PairT x y       -> tupToStr [show x, show y]
      Tup3T x y z     -> tupToStr [show x, show y, show z]
      Tup4T w x y z   -> tupToStr [show w, show x, show y, show z]
      Tup5T v w x y z -> tupToStr [show v, show w, show x, show y, show z]
    where
      tupToStr xs = concat ["Tuple(", intercalate ", " xs, ")"]

--------------------------------------------------------------------------------

-- | Get all 'Variable's in a 'Tag'.
varsTag :: Tag a -> [Variable]
varsTag tg = case tg of
  NullT           -> []
  TrueT           -> []
  FalseT          -> []
  BoolT x y       -> [x, y]
  RealT eq        -> varsEq eq
  CondT x y z     -> nub $ concat [varsTag x, varsTag y, varsTag z]
  PairT x y       -> nub $ varsTag x ++ varsTag y
  Tup3T x y z     -> nub $ concat [varsTag x, varsTag y, varsTag z]
  Tup4T w x y z   -> nub $ concat [varsTag w, varsTag x, varsTag y, varsTag z]
  Tup5T v w x y z -> nub $ concat [varsTag v, varsTag w, varsTag x, varsTag y, varsTag z]

-- | Get all 'Variable's in an 'Equation'.
varsEq :: Equation -> [Variable]
varsEq eq = nub $ concatMap vars eq

--------------------------------------------------------------------------------

-- | Apply the given function to each variable in the 'Tag'.
applyFtoTag :: (Variable -> Variable) -> Tag a -> Tag a
applyFtoTag f tg = case tg of
  NullT           -> NullT
  TrueT           -> TrueT
  FalseT          -> FalseT
  BoolT x y       -> BoolT (f x) (f y)
  RealT eq        -> RealT (eqF eq f)
  PairT x y       -> PairT (fts x) (fts y)
  Tup3T x y z     -> Tup3T (fts x) (fts y) (fts z)
  Tup4T x y z w   -> Tup4T (fts x) (fts y) (fts z) (fts w)
  Tup5T x y z w u -> Tup5T (fts x) (fts y) (fts z) (fts w) (fts u)
  CondT x y z     -> CondT (fts x) (fts y) (fts z)
  where
    fts :: Tag a -> Tag a
    fts s = applyFtoTag f s
    eqF eq fun = [Term c [fun x | x <- v] | (Term c v) <- eq]

--------------------------------------------------------------------------------

-- | Flag that tracks the state of compilation of a variable.
--
-- Uncompiled: Variable is still single-railed and positive or negative
-- Compiled: Varible is dual-railed but it's final ODE is not yet split between
--           its positive and negative parts.
-- Finalized: Variable is dual-railed and it's ODEs are already finalized
data CompileFlag = Uncompiled | Compiled | Finalized deriving (Eq)

-- | The provenance or origin of an ODE are a set of tags denoting where a
--   particular reaction came from.
type Origin = S.Set String

origin :: String -> Origin
origin = S.singleton

-- Every variable has an ODE, which includes it's derivative equation,
-- it's initial condition, and it's current state of compilation
data ODE = ODE { equation :: Equation
               , initial  :: Double
               , flag     :: CompileFlag
               , source   :: Maybe Origin
               }

-- | A \"typed\" signal represented by a system of ODEs.
data Signal a = Sg { tag    :: Tag a -- ^ /Type/ of the SF
                   , system :: [ODE] -- ^ System of equations
                   }

-- | A signal function, or SF transformer.
newtype SF a b = SF { runSF :: Signal a -> Signal b }

fromSource :: String -> SF a a
fromSource src = SF $ \(Sg tg sys) ->
  Sg tg (map (\(ODE eqn i flg _) -> ODE eqn i flg (Just $ origin src)) sys)

--------------------------------------------------------------------------------

-- | Identity: emit the input signal unaltered.
idSF :: SF a a
idSF = SF id

-- | Lift a pure tg function to an SF.
--
-- You can regard the tag function as a "reindexing" the structure of the
-- tag object without actually modifying the underlying system of equations
-- or initial conditions. The resulting SF is really an "adapter" that makes two
-- other SFs compatible.
arrTag :: (Tag a -> Tag b) -> SF a b
arrTag f = SF $ \(Sg tg sys) -> Sg (f tg) sys

-- | Extract the first element of a signal.
proj1SF :: SF (a, b) a
proj1SF = arrTag $ \(PairT x _) -> x

-- | Extract the second element of a signal.
proj2SF :: SF (a, b) b
proj2SF = arrTag $ \(PairT _ y) -> y

--------------------------------------------------------------------------------

-- | Create an SF that ignores its input, producing the constant signal
-- s regardless.
constSF :: Signal a -> SF b a
constSF s = SF (const s)

-- | Duplicate: emit a pair of the input signal.
dupSF :: SF a (a, a)
dupSF = arrTag $ \x -> PairT x x

--------------------------------------------------------------------------------

-- | Convert a 'Tup3T' signal into a 'PairT' signal.
tup3ToPairSF :: SF (a, b, c) (a, (b, c))
tup3ToPairSF = arrTag $ \(Tup3T x y z) -> PairT x (PairT y z)

-- | Convert a 'Tup4T' signal into a 'PairT' signal.
tup4ToPairSF :: SF (a, b, c, d) (a, (b, c, d))
tup4ToPairSF = arrTag $ \(Tup4T x y z w) -> PairT x (Tup3T y z w)

-- | Convert a 'Tup5T' signal into a 'PairT' signal.
tup5ToPairSF :: SF (a, b, c, d, e) (a, (b, c, d, e))
tup5ToPairSF = arrTag $ \(Tup5T x y z w u) -> PairT x (Tup4T y z w u)

--------------------------------------------------------------------------------

-- | Rewire multiple variables.
rewireVars :: Signal a -> [(Variable, Variable)] -> Signal a
rewireVars s = foldr rewireVar s . modTups
  where
    modTup (x,y) (a,b) = if a == b then (x,y) else (shift x a b, shift y a b)
    modTups tups = foldr (\x acc -> foldr (flip modTup) x acc:acc) [] tups
    shift x a b | x == b = if a > b then a-1 else a
                | x > b = x-1
                | otherwise = x

-- | Rewire 'Variable' d to 'Variable' c.
rewireVar :: (Variable, Variable) ->  Signal a -> Signal a
rewireVar (c, d) sg@(Sg tg sys) | c == d = sg
                                | otherwise = dropVar (Sg (f tg) sys') d
  where
    sys' = [ODE [Term con [replaceVar c d v | v <- var] | (Term con var) <- eq] ic ct src | ODE eq ic ct src <- sys]

    replaceVar :: Variable -> Variable -> Variable -> Variable
    replaceVar cv dv x = if x == dv then cv else x
    f = applyFtoTag (replaceVar c d)

--------------------------------------------------------------------------------

-- | A reaction in a SF.
data Reaction = Rn { reactants :: [Variable] -- ^ Reactants
                   , products  :: [Variable] -- ^ Products
                   , constant  :: Double -- ^ Rate constant
                   }
                deriving Eq

-- | Reaction displayed similarly to: @x1 + x2 --{5}-> x3 + x3@
instance Show Reaction where
  show (Rn r p k) = toSt r ++ " ->[" ++ show k ++ "] " ++ toSt p
    where toSt = intercalate " + " . map (("x" ++) . show)

instance Ord Reaction where
  compare (Rn r1 _ k1) (Rn r2 _ k2) = compare (sort r1,k1) (sort r2,k2)

-- | Update a system of 'Equation's to include the effects of a list of
-- 'Reaction's.
applyRns :: [Reaction] -> [ODE] -> [ODE]
applyRns rns sys = map (normalizeODE . extendEq) [0..length sys - 1]
  where extendEq n = ODE (equation (sys !! n) ++ map (toTerm n) rns) (initial (sys !! n)) (flag (sys !! n)) (source (sys !! n))
        toTerm x (Rn r p k) = Term (k * fromIntegral (count x p - count x r)) r
        count x = length . filter (==x)

-- | Simplify the 'Equation' so that all terms are combined and sorted and all
-- zero terms are dropped.
normalizeEq :: Equation -> Equation
normalizeEq = filter notZero . combineTerms . sort . map sortVars
  where notZero (Term c _) = c /= 0
        sortVars (Term c vs) = Term c (sort vs)
        combineTerms [] = []
        combineTerms [t] = [t]
        combineTerms ( t1@(Term c1 v1) : t2@(Term c2 v2) : ts)
          | v1 == v2  = combineTerms $ Term (c1 + c2) v1 : ts
          | otherwise = t1 : combineTerms (t2 : ts)

-- | Simplify the ODE's equation
normalizeODE :: ODE -> ODE
normalizeODE (ODE eq ic ct src) = ODE (normalizeEq eq) ic ct src

--------------------------------------------------------------------------------

-- | Find all the variables that have a net effect on the tg
--
-- This is accomplished by recursively iterating over the ODEs adding variables
-- if they have a net effect on the variables in the tg
depsTag :: Tag a -> [ODE] -> [Variable]
depsTag tg sys =
  let n = length sys
      deps = map (\i -> varsEq (equation (sys !! i))) [0..n-1]
      f vs = let vs' = sort . nub $ concat (vs : [ deps !! v | v <- vs])
             in if vs == vs' then vs else f vs'
  in f (sort . nub $ varsTag tg)

-- | Remove unnecessary variables in the signal.
reduceSF :: [Variable] -> SF a a
reduceSF reserved = SF $ \s@(Sg tg sys) ->
  let used = depsTag tg sys ++ reserved
      notUsed = filter (not . (`elem` used)) [0..length sys-1]
  in dropVars s notUsed

-- | Remove all of the given variables from the signal.
--
-- WARNING: This function simply removes the relevant equations and
--          initial conditions from the signal. It is possible that
--          the variables are used in other equations.
dropVars :: Signal a -> [Variable] -> Signal a
dropVars s = foldl dropVar s . reverse . sort

-- | Remove the given variable from the signal.
--
-- WARNING: This function simply removes the relevant equations and
--          initial conditions from the signal. It is possible that
--          the variables are used in other equations.
dropVar :: Signal a -> Variable -> Signal a
dropVar s v =
  let Sg tg sys = shiftSg v (-1) s
      sys' = dropIndex v sys
  in Sg tg sys'

-- | Remove the element at index n from a list.
dropIndex :: Int -> [a] -> [a]
dropIndex n xs = lft ++ rgt
  where (lft, _:rgt) = splitAt n xs

-- | Shift variables whose index is @>= n@ by @+k@ in the tg and system of
-- equations of the signal.
--
-- For example, @shiftSg 4 3 s@  will shift signal @s@ by modifying all
-- variable indices starting at @4..@ by @+3@.
--
-- __WARNING__: This does not verify the number of equations is the same
-- as the indices of the variables! Use with caution!
shiftSg :: Int -> Int -> Signal a -> Signal a
shiftSg n k (Sg tg sys) = Sg (shiftTag n k tg) (shiftSys n k sys)

-- | Shift variables whose index is @>= n@ by @+k@.
shiftVar :: Int -> Int -> Variable -> Variable
shiftVar n k var = if var < n then var else var + k

-- | Shift variables whose index is @>= n@ by @+k@.
shiftTerm :: Int -> Int -> Term -> Term
shiftTerm n k (Term c v) = Term c (map (shiftVar n k) v)

-- | Shift variables whose index is @>= n@ by @+k@.
shiftEq :: Int -> Int -> Equation -> Equation
shiftEq n k = map (shiftTerm n k)

shiftODE :: Int -> Int -> ODE -> ODE
shiftODE n k (ODE eq ic ct src) = ODE (shiftEq n k eq) ic ct src

-- | Shift variables whose index is @>= n@ by @+k@.
shiftSys :: Int -> Int -> [ODE] -> [ODE]
shiftSys n k = map (shiftODE n k)

-- | Shift variables whose index is @>= n@ by @+k@.
shiftTag :: Int -> Int -> Tag a -> Tag a
shiftTag n k tg = case tg of
  NullT           -> NullT
  TrueT           -> TrueT
  FalseT          -> FalseT
  BoolT x y       -> BoolT (shiftVar n k x) (shiftVar n k y)
  RealT eq        -> RealT (shiftEq n k eq)
  PairT x y       -> PairT (f x) (f y)
  Tup3T x y z     -> Tup3T (f x) (f y) (f z)
  Tup4T x y z w   -> Tup4T (f x) (f y) (f z) (f w)
  Tup5T x y z w u -> Tup5T (f x) (f y) (f z) (f w) (f u)
  CondT x y z     -> CondT (f x) (f y) (f z)
  where
    f :: Tag a -> Tag a
    f = shiftTag n k
