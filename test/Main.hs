module Main where

import Test.Tasty
import Test.Tasty.HUnit

import Bio.Reactamole

import Bool

main :: IO ()
main = defaultMain boolTests
