# Glauber algorithm for the 1D Ising chain

A small [Julia](https://julialang.org/) project that simulates Glauber
single-spin-flip dynamics on a 1D Ising chain with periodic boundary
conditions, optionally in the presence of a uniform external field `h`.
Used as a test for see how Claude code is working.

Authored by Leandro M. Chinellato <chinellato.leandro@gmail.com>.

## Layout

```
src/Glauber.jl        # module: energy, magnetization, glauber_alg, ...
scripts/run_sweeps.jl # driver that regenerates the four figures
test/runtests.jl      # correctness tests
Project.toml          # dependency declaration
```

## Reproducing the figures

The project uses a `Project.toml` so the environment is reproducible without
relying on your global Julia setup.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/run_sweeps.jl
```

This regenerates `energy.png`, `magnetization.png`, `energy_with_h.png`, and
`magnetization_with_h.png` at the repository root. A fixed random seed is set
in the script so the output is deterministic.

## Running the tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Using the module interactively

```julia
julia> using Pkg; Pkg.activate(".")
julia> include("src/Glauber.jl"); using .Glauber

julia> σ = generate_random_conf(100);
julia> Es, Ms = glauber_alg(1.0, 1.5, σ; Nsteps = 500);        # h = 0
julia> Es, Ms = glauber_alg(1.0, 1.5, σ; h = 0.3, Nsteps = 500);
```
