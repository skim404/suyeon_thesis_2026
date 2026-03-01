include("spatial_model.jl")
using Plots

## Setup
seed = 3

grid = GridSpace(5, 5)
initialize_grid!(grid,
    population_fractions = [0.5, 0.5],
    population_types = [
        [0, 100, 0, 0, 0, 0, 0, 0],
        [0, 0, 100, 0, 0, 0, 0, 0]
    ],
    seed = seed
)
update_totals!(grid)

params = ModelParameters(
    1.0,          ## λ
    0.925,        ## λP
    0.95,         ## λR
    0.87875,      ## λRP
    0.8348125,    ## λRRP
    0.005, 0.2,   ## s_ss, β_ss
    0.0055, 0.25, ## s1, β1
    0.005,  0.2,  ## s2, β2
    1.0,          ## α
    1.0,          ## dose
    1000.0,      ## K
    1.0e-4,       ## migration_rate
    1.0e-4,       ## extinction_rate
    100.0,        ## dt (Gillespie window per step)
    300.0,        ## AB_period
    330.0,        ## no_AB_period
    true,         ## out_of_phase?
    0.5           ## phase_fraction
)

sim = GridDynamics(grid, params; seed = seed)

## Run
println("Running spatial simulation...")
@time run_simulation!(sim; n_steps = 50)

## Summary
ts         = sim.time_series
cell_types = [:nS0, :nSR1, :nSR2, :nSS, :nR0, :nRS, :nRR1, :nRR2]

println("\nFinal average density per patch:")
for col in cell_types
    val = getfield(ts, col)[end]
    println("  $col : $(round(val, digits=2))")
end

## Plot and save

p1 = Plots.plot(get_color_matrix(sim.grid),
                axis=false, title="Spatial distribution")

ys   = hcat(ts.nS0, ts.nSR1, ts.nSR2, ts.nSS, ts.nR0, ts.nRS, ts.nRR1, ts.nRR2)
labs = ["nS0" "nSR1" "nSR2" "nSS" "nR0" "nRS" "nRR1" "nRR2"]
cols = [:gray :blue :green :lightblue :orange :red :purple :brown]

p2 = Plots.plot(ts.times, ys,
                label=labs, color=cols,
                xlabel="Time", ylabel="Average cell count per patch",
                title="")

p = Plots.plot(p1, p2, layout=(1, 2), size=(1100, 450))

savefig(p, joinpath(@__DIR__, "timeseries.png"))
println("\nPlot saved to timeseries.png")
display(p)