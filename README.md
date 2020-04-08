# WCA1

Code for the paper [_A bio-inspired geometric model for sound reconstruction_](https://arxiv.org/abs/2004.02450) by Ugo Boscain, Dario Prandi, Ludovic Sacchelli, Giuseppina Turco.

## How to reproduce the experiments of the paper

From inside the `experiments` folder, run

```sh
julia run_experiments.jl
```

## Usage

After cloning the repository, import the package via

```julia
import Pkg
Pkg.activate(path_to_WCA1_folder)

using WCA1
```

Look at the `notebooks/` folder for an usage example.