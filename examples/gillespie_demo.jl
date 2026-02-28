using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using PlasmidModels
using .PlasmidModels.GillespieModel

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
    40.0,       # AB_period
    350.0       # no_AB_period
)

x0 = Int32[0, 1000, 0, 0, 1000, 0, 0, 0]  # [S0, SR, SR2, SS, R0, RS, RR, RR2]

df = run_gillespie_patch(x0, params; dt=1.0, t_end=500.0, seed=3,
                         output_dir="output/gillespie", save_csv=true)

println("Gillespie demo done. Rows: ", nrow(df))
