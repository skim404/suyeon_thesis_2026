## ode_standalone.jl
using DifferentialEquations, CSV, DataFrames, Plots, Statistics

struct SimulationParams{T<:Real}
    λ::T
    α::T
    β1::T      # Conjugation rate for SR (resistant plasmid 1)
    s1::T      # Segregation rate for SR (and SS, the sensitive plasmid)
    x1::T      # Base plasmid cost for SR and SS
    β2::T      # Conjugation rate for SR2 (resistant plasmid 2)
    s2::T      # Segregation rate for SR2
    x2::T      # Base plasmid cost for SR2
    cR::T      # Chromosomal resistance cost
    x_slope::T # Slope for conjugation-dependent plasmid cost
end

plasmid_cost(β::T, x_base::T, x_slope::T) where {T<:Real} = x_base + x_slope * β

function vector_field!(du, u, p::SimulationParams, t)
    # u = [S0, SR, SR2, SS, R0, RS, RR, RR2, A]
    nS0, nSR, nSR2, nSS, nR0, nRS, nRR, nRR2, A = u

    λ, α, β1, s1, x1_base, β2, s2, x2_base, cR, x_slope =
        p.λ, p.α, p.β1, p.s1, p.x1, p.β2, p.s2, p.x2, p.cR, p.x_slope

    x1 = plasmid_cost(β1, x1_base, x_slope)  # for SR and SS
    x2 = plasmid_cost(β2, x2_base, x_slope)  # for SR2

    N = nS0 + nSR + nSR2 + nSS + nR0 + nRS + nRR + nRR2
    αN = α * N

    one_minus_cR = 1 - cR

    one_minus_x1 = 1 - x1
    λSS  = one_minus_x1 * λ
    λRR  = one_minus_cR * one_minus_x1 * λ    # growth for chrom. resistant + plasmid1/sensitive-plasmid background

    one_minus_x2 = 1 - x2
    λRR2 = one_minus_cR * one_minus_x2 * λ    # growth for chrom. resistant + plasmid2 background

    λR = one_minus_cR * λ                      # chrom. resistant, plasmid-free

    # total conjugation pressure into plasmid-free recipients
    β_sum = β1 * (nSR + nRR) + β2 * (nSR2 + nRR2) + β1 * (nSS + nRS)

    # S0
    du[1] = nS0 * (λ - β_sum - αN - A) +
            s1 * (λSS * nSS + λRR * nSR) +
            s2 * (λRR2 * nSR2)

    # SR
    du[2] = nSR * (λRR * (1 - s1) - αN) +
            β1 * (nSR + nRR) * nS0

    # SR2
    du[3] = nSR2 * (λRR2 * (1 - s2) - αN) +
            β2 * (nSR2 + nRR2) * nS0

    # SS (sensitive plasmid)
    du[4] = nSS * (λSS * (1 - s1) - αN - A) +
            β1 * (nSS + nRS) * nS0

    # R0
    du[5] = nR0 * (λR - β_sum - αN) +
            s1 * (λRR * nRS + λRR * nRR) +
            s2 * (λRR2 * nRR2)

    # RS
    du[6] = nRS * (λRR * (1 - s1) - αN) +
            β1 * (nSS + nRS) * nR0

    # RR
    du[7] = nRR * (λRR * (1 - s1) - αN) +
            β1 * (nSR + nRR) * nR0

    # RR2
    du[8] = nRR2 * (λRR2 * (1 - s2) - αN) +
            β2 * (nSR2 + nRR2) * nR0

    # A is externally controlled
    du[9] = 0.0
    return nothing
end

function run_simulation(tspan, u0, no_AB_period, AB_period, params::SimulationParams)
    all_ts = Float64[]
    all_ys = Vector{Vector{Float64}}()

    current_t = float(tspan[1])
    current_y = Float64.(u0)

    while current_t < tspan[2]
        # no antibiotic
        end_t = min(current_t + no_AB_period, tspan[2])
        current_y[9] = 0.0
        prob = ODEProblem(vector_field!, current_y, (current_t, end_t), params)
        sol = solve(prob, Tsit5(); adaptive=false, dt=0.2)
        append!(all_ts, sol.t)
        append!(all_ys, sol.u)
        current_y .= sol.u[end]
        current_t = sol.t[end]
        current_t >= tspan[2] && break

        # antibiotic present
        end_t = min(current_t + AB_period, tspan[2])
        current_y[9] = 1.0
        prob = ODEProblem(vector_field!, current_y, (current_t, end_t), params)
        sol = solve(prob, Tsit5(); adaptive=false, dt=0.2)
        append!(all_ts, sol.t)
        append!(all_ys, sol.u)
        current_y .= sol.u[end]
        current_t = sol.t[end]
    end

    col_names = [:time, :nS0, :nSR, :nSR2, :nSS, :nR0, :nRS, :nRR, :nRR2, :A]
    data = hcat(all_ts, reduce(hcat, all_ys)')
    return DataFrame(data, col_names)
end

function single_simulation(u0, tspan, params::SimulationParams, no_AB, AB;
                          output_dir="output_ode", save_csv=true, save_plots=true)
    isdir(output_dir) || mkpath(output_dir)

    results = run_simulation(tspan, u0, no_AB, AB, params)

    if save_csv
        CSV.write(joinpath(output_dir, "simulation.csv"), results)
    end

    if save_plots
        p_full = plot(results.time,
            [results.nS0 results.nSR results.nSR2 results.nSS results.nR0 results.nRS results.nRR results.nRR2],
            xlabel="Time", ylabel="Population density",
            title="ODE simulation",
            label=["S0" "SR" "SR2" "SS" "R0" "RS" "RR" "RR2"],
            linewidth=2, legend=:outerright)

        p_log = plot(results.time,
            [results.nS0 results.nSR results.nSR2 results.nSS results.nR0 results.nRS results.nRR results.nRR2],
            xlabel="Time", ylabel="Population density (log10)",
            title="ODE simulation (log scale)",
            label=["S0" "SR" "SR2" "SS" "R0" "RS" "RR" "RR2"],
            yscale=:log10, ylims=(1e-12, 1.0),
            linewidth=2, legend=:outerright)

        savefig(p_full, joinpath(output_dir, "full.png"))
        savefig(p_log,  joinpath(output_dir, "full_log.png"))
    end

    return results
end

# =========================
# MAIN
# =========================
println("\n" * "="^70)
println("ODE: SINGLE RUN")
println("="^70)

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

@time results = single_simulation(u0, tspan, params, no_AB_period, AB_period;
                                 output_dir="output_ode", save_csv=true, save_plots=true)

println("Done. Rows: ", nrow(results))
