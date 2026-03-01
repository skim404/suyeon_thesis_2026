include("ode_model.jl")
## Parameters

## nSS, nSR1, nSR2, nSS, nR0, nRS, nRR1, nRR2
u0    = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0]
tspan = (0.0, 3000.0)
args  = (
    1.0,   ## λ  : base growth rate
    1.0,   ## α  : death rate
    0.2,   ## β1 : conjugation rate, plasmid 1
    0.2,   ## β2 : conjugation rate, plasmid 2
    0.005, ## s1 : segregation rate, plasmid 1
    0.005, ## s2 : segregation rate, plasmid 2
    0.075, ## x  : plasmid carriage cost
    0.05   ## cR : resistance cost
)
no_AB = 150.0
AB    = 150.0

## Run
println("Running ODE simulation...")
@time results = run_simulation(tspan, u0, no_AB, AB, args)

## Summary
cell_types = [:nS0, :nSR1, :nSR2, :nSS, :nR0, :nRS, :nRR1, :nRR2]
final      = results[end, :]

println("\nFinal densities:")
for col in cell_types
    println("  $col : $(round(final[col], digits=6))")
end

## Plot and save
ys   = Matrix(results[:, [:nS0, :nSR1, :nSR2, :nSS, :nR0, :nRS, :nRR1, :nRR2]])
labs = ["nS0" "nSR1" "nSR2" "nSS" "nR0" "nRS" "nRR1" "nRR2"]
cols = [:gray :blue :green :lightblue :orange :red :purple :brown]

p = Plots.plot(results.time, ys,
               label=labs, color=cols,
               xlabel="Time", ylabel="Population density",
               title="",
               ylims=(0, 1.0))

savefig(p, joinpath(@__DIR__, "timeseries.png"))
println("\nPlot saved to timeseries.png")
display(p)