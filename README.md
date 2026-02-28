# Evolution of Plasmid-Carried Antibiotic Resistance Genes in Fluctuating Environments

1. **ODE model** 
2. **Spatial model** 
3. **Gillespie (single patch)**

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

