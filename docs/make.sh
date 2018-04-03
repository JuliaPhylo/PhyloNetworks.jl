#!/bin/bash
# Idea of this script taken from: http://steven.casagrande.io/articles/travis-ci-and-if-statements/

set -ev

if [ "$TRAVIS_OS_NAME" == "linux" ]; then
    julia -e 'Pkg.clone("https://github.com/pbastide/Documenter.jl")';
    #julia -e 'Pkg.checkout("Documenter", "update_documenter")';
    julia -e 'Pkg.clone("https://github.com/cecileane/PhyloPlots.jl")';
    #julia -e 'Pkg.add("Cairo")';
    #julia -e 'Pkg.add("Fontconfig")';
    #rm $HOME/.julia/lib/v$TRAVIS_JULIA_VERSION/Compose.ji;
    julia -e 'Pkg.add("Weave")';
    julia -e 'cd(Pkg.dir("PhyloNetworks")); include(joinpath("docs", "make.jl"))';
fi

exit 0;
