using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using PlasmidModels
using .PlasmidModels.ODEModel

λ = 1.0
α = 1.0
β1 = 0.12
s1 = 0.05
β2 = 0.115
s2 = 0.05
x_base = 0.075
x_slope = 0.1
cR = 0.05

params = SimulationParams(λ, α, β1, s1, x_base, β2, s2, x_base, cR, x_slope)

u0 = [0.0, 0.9999, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # [S0, SR, SR2, SS, R0, RS, RR, RR2, A]
tspan = (0.0, 5000.0)

no_AB_period = 350.0
AB_period = 40.0

results = single_simulation(u0, tspan, params, no_AB_period, AB_period;
                            output_dir="output/ode",
                            save_csv=true,
                            save_plots=true)

println("ODE demo done. Rows: ", nrow(results))
