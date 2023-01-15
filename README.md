# Pre-compiling

in REPL:

```julia
using PackageCompiler
create_sysimage(sysimage_path="sysimage_ardca.dylib", precompile_execution_file="scripts/run_gdca_precompile.jl")
```

To see the difference, compare:
```julia
julia --trace-compile=stderr --sysimage sysimage_ardca.dylib scripts/run_gdca_precompile.jl
julia --trace-compile=stderr scripts/run_gdca_precompile.jl
```

n.b. here we omit the packages argument to create_sysimage.
This ensures that all packages in the project are put in the sysimage.

N.B. the sysimage can itself be redistributed - as long as packages don't
rely on global paths to files on the local filesystem:
https://stackoverflow.com/questions/69556287/moving-a-precompiled-sysimage-to-a-new-machine

TODO think about making an App. The only difficulty is we need to probably
move from using ArgParse to parsing cmd line args via ARGS
so we should write some new entrypoints that handle this.
https://julialang.github.io/PackageCompiler.jl/dev/apps.html#Creating-an-app
https://github.com/JuliaLang/PackageCompiler.jl/tree/master/examples/MyApp

Ref:
https://www.youtube.com/watch?v=d7avhSuK2NA


# ArDCA

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pagnani.github.io/ArDCA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pagnani.github.io/ArDCA.jl/dev)
[![Build Status](https://github.com/pagnani/ArDCA/workflows/CI/badge.svg)](https://github.com/pagnani/ArDCA/actions)
[![Coverage](https://codecov.io/gh/pagnani/ArDCA/branch/master/graph/badge.svg)](https://codecov.io/gh/pagnani/ArDCA)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


Autoregressive protein model learning through generalized logistic regression in Julia.
## Overview

The authors of this code are Jeanne Trinquier, Guido Uguzzoni, Andrea Pagnani, Francesco Zamponi, and Martin Weigt.

See also [this Wikipedia article](https://en.wikipedia.org/wiki/Direct_coupling_analysis) article for a general overview of the Direct Coupling Analysis technique. 

The code is written in [Julia](https://www.julialang.org/).

## Install

This is a registered package: to install enter `]` in the repl and

```
pkg> add ArDCA 
```
## Notebooks

There are two `jupyter` notebooks (Python, and Julia) to help using the Package.

The [tutorial.ipynb](julia-notebook/tutorial.ipynb) is for the julia version.
The [arDCA_sklearn.ipynb](python-notebook/arDCA_sklearn.ipynb) is for the python version.

## Data 

Data for five protein families (PF00014,PF00072, PF00076,PF00595,PF13354) are contained in the companion
[ArDCAData](https://github.com/pagnani/ArDCAData) package.

For didactic reasons we include locally in the `data` folder, the PF00014 dataset.

## Requirements

The minimal Julia version to run this code is 1.5. To run it in parallel 
using Julia multicore infrastructure, start julia with

```
$> julia -t numcores # ncores can be as large as your available number of threads
```

## Documentation

[Stable version](https://pagnani.github.io/ArDCA.jl/stable)

[Development version](https://pagnani.github.io/ArDCA.jl/dev)

## License

This project is covered under the MIT License.


## Example of loading and saving model

using ArDCA
arnet,arvar = ardca("/Users/alex/proteins/aflatent/data/family/cm/cm_match_msa.fa", output_file="inttest")
loglikelihood("TSENPLLALREKISALDEKLLALLAERRELAVEVGKAKLLSHRPVRDIDRERDLLERLITLGKAHHLDAHYITRLFQLIIEDSVLTQQALLQQH", arnet)
arnetreload = load_arnet("inttest")
loglikelihood("TSENPLLALREKISALDEKLLALLAERRELAVEVGKAKLLSHRPVRDIDRERDLLERLITLGKAHHLDAHYITRLFQLIIEDSVLTQQALLQQH", arnetreload)