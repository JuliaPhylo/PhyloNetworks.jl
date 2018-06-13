#!/bin/bash
# Idea of this script taken from: http://steven.casagrande.io/articles/travis-ci-and-if-statements/

set -ev

if [ "$TRAVIS_OS_NAME" == "linux" ]; then
    julia -e 'Pkg.clone("https://github.com/cecileane/PhyloPlots.jl")';
    julia -e 'Pkg.add("Documenter")';
    julia -e 'cd(Pkg.dir("PhyloNetworks")); include(joinpath("docs", "make.jl"))';
fi

exit 0;
