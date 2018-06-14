#!/bin/bash
# Idea of this script taken from: http://steven.casagrande.io/articles/travis-ci-and-if-statements/

set -ev

if [ "$TRAVIS_OS_NAME" == "linux" ]; then
    ## Install R
    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9;
    sudo add-apt-repository -y "deb http://cran.rstudio.com/bin/linux/ubuntu $(lsb_release -s -c)/";
    sudo apt-get update -qq -y;
    sudo apt-get install git r-base r-base-dev r-recommended -y;
    ## Build the doc 
    julia -e 'Pkg.add("PhyloPlots")';
    julia -e 'Pkg.add("Documenter")';
    julia -e 'cd(Pkg.dir("PhyloNetworks")); include(joinpath("docs", "make.jl"))';
fi

exit 0;
