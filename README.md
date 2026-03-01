# Evolution of Plasmid-Carried Antibiotic Resistance Genes in Fluctuating Antibiotic Environments

This project assesses whether plasmid-carried antibiotic resistance genes can be favoured over chromosome-carried ones despite the inherent instability of plasmids in fluctuating antibiotic environments. The ordinary differential equation (ODE) model is adapted from Lehtinen et al. (2019), which was then extended into a Gillespie algorithm model and a spatial model. All three models simulate bacterial populations comprised of 6 (+2, if an additional resistant plasmid is considered), each defined by chromosomal resistance status, plasmid carriage, and plasmid resistance type. These populations are subjected to a fluctuating antibiotic environment whose temporal structure is defined by the two parameters Ton and Toff.

## Models

### [ODE](models/ODE/)
A deterministic, continuous model comprised of ODEs that each describe the change in population density of bacterial cell types. Adapted from Lehtinen et al. (2019). Refer to Section 2.1 of the thesis for further details on the model. This model was used to obtain the results for Section 3.1, 3.3, 3.4, 3.5, 3.6, and 3.7.

### [Gillespie](models/Gillespie/)
Stochastic, discrete model using the Gillespie algorithm (Gillespie, 1977). This model captures demographic stochasticity that can lead to extinction of cell types. Refer to Section 2.2 of the thesis for further details on the model. This model was used to obtain the results for Section 3.2, and it is also incorporated into the spatial model. 

### [Spatial](models/Spatial/)
Spatial extension of the Gillespie model, where each patch on a 2D grid runs an independent Gillespie simulation to describe within-group dynamics, coupled through migration and local extinction that gives rise to between-group dynamics. Refer to Section 2.3 for further details on the model. This model was used to obtain the results for Section 3.8.

## Cell types

All models use the same 8 genotypes:

| Cell Type | Chromosome | Plasmid        |
|----------|-----------|----------------|
| `S0`    | Sensitive  | None           |
| `SR1`   | Sensitive  | Resistant P1   |
| `SR2`   | Sensitive  | Resistant P2   |
| `SS`    | Sensitive  | Sensitive      |
| `R0`    | Resistant  | None           |
| `RS`    | Resistant  | Sensitive      |
| `RR1`   | Resistant  | Resistant P1   |
| `RR2`   | Resistant  | Resistant P2   |

## Prerequisites

- Julia ≥ 1.9

Install dependencies from the repo root:

```julia
using Pkg
Pkg.instantiate()
```

## Running the Demo

```
# ODE model
julia models/ODE/demo.jl

# Gillespie model
julia models/Gillespie/demo.jl

# Spatial model (multithreading recommended)
julia --threads auto models/Spatial/demo.jl
```
