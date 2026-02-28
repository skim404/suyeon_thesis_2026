using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using PlasmidModels
using .PlasmidModels.SpatialModel

grid = GridSpace(20, 20)
initialize_grid!(grid;
    population_fractions = [0.5, 0.25, 0.25],
    population_types = [
        [3000, 3000, 0,    0,    0,    0,    0, 0],   # S0 + SR
        [3000, 0,    3000, 0,    0,    0,    0, 0],   # S0 + SR2
        [0,    0,    0,    0, 3000, 3000,    0, 0]    # R0 + RS
    ],
    seed=3
)

params = ModelParameters(
    1.0,        # λ
    0.925,      # λSS
    0.95,       # λR
    0.87875,    # λRR
    0.8348125,  # λRR2
    0.005, 0.2, # SS_s, SS_β
    0.005, 0.2, # s1, β1
    0.032, 0.14,# s2, β2
    1.0,        # α
    1.0,        # dose
    1_000_000.0,# K
    0.001,      # migration_rate
    0.0,        # extinction_rate
    1.0,        # dt
    40.0,       # AB_period
    350.0,      # no_AB_period
    true,       # out_of_phase
    0.5         # phase_fraction
)

sim = GridDynamics(grid, params; seed=3)
run_simulation!(sim; n_steps=200, display_interval=50, output_dir="output/spatial_frames")

println("Spatial demo done.")
