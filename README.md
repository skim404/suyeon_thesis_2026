# PlasmidModels (Julia)

A Julia repo packaging three bacterial population models:

1. **ODE model** (well-mixed; two resistant plasmids + one sensitive plasmid; fluctuating antibiotic)
2. **Spatial model** (grid of patches; within-patch SSA + optional migration/extinction + frame output)
3. **Gillespie (single patch)** extracted from the spatial SSA (non-spatial).

Cell-type notation is consistent across models:

| Symbol | Meaning |
|---|---|
| `S0` | Chromosomally sensitive, plasmid-free |
| `SS` | Chromosomally sensitive, **sensitive** plasmid |
| `SR` | Chromosomally sensitive, resistant plasmid 1 |
| `SR2` | Chromosomally sensitive, resistant plasmid 2 |
| `R0` | Chromosomally resistant, plasmid-free |
| `RS` | Chromosomally resistant, sensitive plasmid |
| `RR` | Chromosomally resistant, resistant plasmid 1 |
| `RR2` | Chromosomally resistant, resistant plasmid 2 |

## Install / run demos

Clone + instantiate:
```bash
git clone <YOUR_GITHUB_URL_HERE>
cd PlasmidModels
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Run demos (write outputs into `output/`):

```bash
julia --project=. examples/ode_demo.jl
julia --project=. examples/gillespie_demo.jl
julia --project=. examples/spatial_demo.jl
```

## Using as a library
From a Julia session started in this repo:
```julia
using PlasmidModels
```

Modules:
- `PlasmidModels.ODEModel`
- `PlasmidModels.GillespieModel`
- `PlasmidModels.SpatialModel`

## Outputs
- ODE demo: `output/ode/` (CSVs + PNG plots)
- Gillespie demo: `output/gillespie/` (CSV)
- Spatial demo: `output/spatial_frames/` (PNG frames)

