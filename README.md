# QPNets.jl

[![License](https://img.shields.io/badge/license-MIT-blue)](https://opensource.org/licenses/MIT)

A package for modeling Quadratic Program Networks (QPNets), and computing corresponding equilibrium points. Examples can be found in the examples subdirectory. This package is under active development. Medium and Large-scale problems may be difficult or impossible to solve using the algorithms implemented in this package at this time. Bug reports are welcomed and appreciated!

[The paper introducing QPNets can be found here.](https://arxiv.org/abs/2404.03767) If you use this package in your research, please use the following citation.

```bibtex
@misc{laine2024mathematical,
      title={Mathematical Program Networks}, 
      author={Forrest Laine},
      year={2024},
      eprint={2404.03767},
      archivePrefix={arXiv},
      primaryClass={math.OC}
}
```

# Installation

This package is not yet publically registered. To install, simply clone this repository locally, navigate to the directory, and run the following:

```
julia --project
```
```julia
# hit `]` to enter "pkg"-mode of the REPL
pkg> instantiate
# hit `delete` or `backspace` to exit "pkg"-mode of the REPL
julia> using QPNets
```

# Running examples

Examples from the above-linked paper can be run using the following interface. 

```julia
julia> qpn = setup(:robust_avoid_simple);
julia> ret = solve(qpn);
```

# Constructing new problems

New QPNets can be constructed by following the template provided in the ``` examples ``` directory. Every QPNet should define a ```setup``` function, which defines the decision variables, QP nodes, and edges of the network, and returns a QPNet object. A default initialization can be assigned to the QPNet, in which case an equilibrium can be computed by simply calling ```solve``` on the QPNet. If not specified, the default initialization is the zero vector. A custom initialization can be passed to the solver by calling 
```solve(qpn, init)```.

# Documentation

A formal documentation of this package is in the works.
