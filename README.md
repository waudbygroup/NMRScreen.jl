# NMRScreen

[![Build Status](https://github.com/waudbygroup/NMRScreen.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/waudbygroup/NMRScreen.jl/actions/workflows/CI.yml?query=branch%3Amain)

# NMRScreen.jl

A Julia package for screening NMR experiments.

## Installation

Since this package is not yet registered, the package must be downloaded directly from GitHub.

## Quick Start

1. Navigate to package directory:
```bash
cd /path/to/NMRScreen.jl
```

2. Activate and instantiate environment:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

3. Use the package:
```julia
using NMRScreen
screen("example/experiment.toml")
```

## TODO

- [x] Add heat maps for base R2 and CSPs
- [x] Test with multiple cocktails and different numbers of peaks
- [x] Sort peaks by chemical shift
- [x] left/right buttons
- [x] summary plots (ranked R2)
- [x] save button
- [x] adjust peak positions
- [x] align spectra before calculating CSPs
- [x] fix scaling of spectra maxima
- [x] shift+keys to move points faster
- [ ] click points on spectrum to highlight
- [ ] drag points on spectrum to move
- [x] navigate between multiple peaks for fragment
- [ ] make separate ref/expt directories