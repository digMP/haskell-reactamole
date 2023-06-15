.PHONY : build
build :
	cabal build

.PHONY : clean
clean :
	cabal clean

.PHONY : repl
repl :
	cabal repl

.PHONY : docs
docs :
	cabal haddock
