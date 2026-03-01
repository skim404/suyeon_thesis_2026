include("gillespie_model.jl")
using Plots

params = GillespieParameters(
    1.0,          ## λ    : base growth rate
    0.925,        ## λP   : growth rate with sensitive plasmid
    0.95,         ## λR   : growth rate with chrom. resistance
    0.87875,      ## λRP  : growth rate with resistance + plasmid (sens. chromosome)
    0.87875,    ## λRRP : growth rate with chrom. resistance + plasmid; no double cost
    0.005, 0.2,   ## s_ss, β_ss : sensitive plasmid dynamics
    0.0055, 0.25, ## s1, β1     : resistant plasmid 1
    0.005,  0.2,  ## s2, β2     : resistant plasmid 2
    1.0,          ## α    : death rate
    1.0,          ## dose : antibiotic dose
    10000.0,      ## K    : carrying capacity
    150.0,        ## AB_period
    150.0         ## no_AB_period
)

## Initial state [nS0, nSR1, nSR2, nSS, nR0, nRS, nRR1, nRR2]
x0 = Int32[0, 1000, 1000, 0, 0, 0, 0, 0]

## Run
println("Running Gillespie...")
@time results = run_gillespie(x0, params, 3000.0;
                              record_interval = 50.0,
                              phase = Int8(0),
                              seed  = 42)

## Summary
cell_types = [:nS0, :nSR1, :nSR2, :nSS, :nR0, :nRS, :nRR1, :nRR2]
final      = results[end, :]

println("\nFinal cell counts:")
for col in cell_types
    println("  $col : $(round(Int, final[col]))")
end

## Plot
ys   = Matrix(results[:, [:nS0, :nSR1, :nSR2, :nSS, :nR0, :nRS, :nRR1, :nRR2]])
labs = ["nS0" "nSR1" "nSR2" "nSS" "nR0" "nRS" "nRR1" "nRR2"]
cols = [:gray :blue :green :lightblue :orange :red :purple :brown]

p = Plots.plot(results.time, ys,
               label=labs, color=cols,
               xlabel="Time", ylabel="Cell count",
               title="")

savefig(p, joinpath(@__DIR__, "timeseries.png"))
println("\nPlot saved to timeseries.png")
display(p)