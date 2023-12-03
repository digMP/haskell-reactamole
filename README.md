# Reactamole: A Functional Reactive Molecular Programming DSL

Reactamole is a domain-specific language for molecular programming that
utilizes reactive functional programming principles. Reactamole observes a
direct correspondence between signal functions in a reactive functional program
(FRP) and chemical reaction networks (CRNs) in a molecular program. This
correspondence allows us to directly translate the core combinators of
(arrowized) FRP to CRNs. Because Reactamole is embedded in the Haskell
programming language, it takes advantage of Haskell's rich type system to
ensure the well-formedness of the resulting chemical reaction networks.

## Installation

Reactamole is implemented as a Haskell eDSL, accessible as a Haskell library.

1.  Download `ghcup` (https://www.haskell.org/ghcup) to install Haskell and its
    toolchain.
2.  `cabal build` to build the project.
3.  `cabal repl` to run Reactamole in GHCi.
4.  `cabal haddock` to build API documentation.

## Running Examples

The `Bio.Reactamole.Examples` module contains a number of example molecular
programs for you to explore. You can use the functions from
`Bio.Reactamole.Export` to export these Reactamole programs as ODEs or
collections of reactions. Here is an example of their usage:

~~~console
$> cabal repl
GHCi, version 9.6.1: https://www.haskell.org/ghc/  :? for help
Loaded GHCi configuration from ...
λ> import Bio.Reactamole
λ> import Bio.Reactamole.Examples
λ> :t srLatch
srLatch :: SF (Bool, Bool) (Bool, Bool)
λ> toIVP srLatch
INPUT:
  Tuple(Bool(x0, x1), Bool(x2, x3))

OUTPUT:
  Tuple(Bool(x4, x5), Bool(x6, x7))

ODEs:
  x4 = 0.0, dx4/dt = [-30.0*x0*x4*x6,+30.0*x1*x5,+90.0*x4^2*x5,-90.0*x4*x5^2,+30.0*x5*x7] [nandSF]
  x5 = 1.0, dx5/dt = [+30.0*x0*x4*x6,-30.0*x1*x5,-90.0*x4^2*x5,+90.0*x4*x5^2,-30.0*x5*x7] [nandSF]
  x6 = 0.0, dx6/dt = [-30.0*x2*x4*x6,+30.0*x3*x7,+30.0*x5*x7,+90.0*x6^2*x7,-90.0*x6*x7^2] [nandSF]
  x7 = 1.0, dx7/dt = [+30.0*x2*x4*x6,-30.0*x3*x7,-30.0*x5*x7,-90.0*x6^2*x7,+90.0*x6*x7^2] [nandSF]

λ> toCRN srLatch
INPUT:
  Tuple(Bool(x0, x1), Bool(x2, x3))

OUTPUT:
  Tuple(Bool(x4, x5), Bool(x6, x7))

REACTIONS:
  x0 + x4 + x6 ->[30.0] x0 + x6 + x5 [nandSF]
  x1 + x5 ->[30.0] x1 + x4 [nandSF]
  x2 + x4 + x6 ->[30.0] x2 + x4 + x7 [nandSF]
  x3 + x7 ->[30.0] x3 + x6 [nandSF]
  x4 + x4 + x5 ->[90.0] x4 + x4 + x4 [nandSF]
  x4 + x5 + x5 ->[90.0] x5 + x5 + x5 [nandSF]
  x5 + x7 ->[30.0] x4 + x6 [nandSF]
  x6 + x6 + x7 ->[90.0] x6 + x6 + x6 [nandSF]
  x6 + x7 + x7 ->[90.0] x7 + x7 + x7 [nandSF]

INITIAL CONDITIONS:
  x4 = 0.0
  x5 = 1.0
  x6 = 0.0
  x7 = 1.0
~~~

## Publication

Titus H. Klinge, James I. Lathrop, Peter-Michael Osera, and Allison Rogers.
Reactamole: Functional Reactive Molecular Programming. 27th International
Conference on DNA Computing and Molecular Programming (DNA '27). September,
2021, Oxford, UK.
