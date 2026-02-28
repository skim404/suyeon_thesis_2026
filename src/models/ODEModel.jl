\
module ODEModel

using DifferentialEquations
using CSV
using DataFrames
using Plots
using Statistics

export SimulationParams, plasmid_cost, vector_field!, run_simulation, single_simulation

struct SimulationParams{T<:Real}
    λ::T
    α::T
    β1::T      # Conjugation rate for SR (resistant plasmid 1)
    s1::T      # Segregation rate for SR
    x1::T      # Base plasmid cost (used for SR and SS)
    β2::T      # Conjugation rate for SR2 (resistant plasmid 2)
    s2::T      # Segregation rate for SR2
    x2::T      # Base plasmid cost (used for SR2)
    cR::T      # Chromosomal resistance cost
    x_slope::T # Slope for conjugation-dependent plasmid cost
end

function plasmid_cost(β::T, x_base::T, x_slope::T) where {T<:Real}
    return x_base + x_slope * β
end

function vector_field!(du, u, p::SimulationParams, t)
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
    λRP1 = one_minus_cR * one_minus_x1 * λ

    one_minus_x2 = 1 - x2
    λRP2 = one_minus_cR * one_minus_x2 * λ

    λR = one_minus_cR * λ

    β_sum = β1 * (nSR + nRR) + β2 * (nSR2 + nRR2) + β1 * (nSS + nRS)

    du[1] = nS0 * (λ - β_sum - αN - A) +
            s1 * (λSS * nSS + λRP1 * nSR) +
            s2 * (λRP2 * nSR2)

    du[2] = nSR * (λRP1 * (1 - s1) - αN) +
            β1 * (nSR + nRR) * nS0

    du[3] = nSR2 * (λRP2 * (1 - s2) - αN) +
            β2 * (nSR2 + nRR2) * nS0

    du[4] = nSS * (λSS * (1 - s1) - αN - A) +
            β1 * (nSS + nRS) * nS0

    du[5] = nR0 * (λR - β_sum - αN) +
            s1 * (λRP1 * nRR + λRP1 * nRS) +
            s2 * (λRP2 * nRR2)

    du[6] = nRS * (λRP1 * (1 - s1) - αN) +
            β1 * (nSS + nRS) * nR0

    du[7] = nRR * (λRP1 * (1 - s1) - αN) +
            β1 * (nSR + nRR) * nR0

    du[8] = nRR2 * (λRP2 * (1 - s2) - αN) +
            β2 * (nSR2 + nRR2) * nR0

    du[9] = 0.0
    return nothing
end

function run_simulation(tspan::Tuple{Real,Real}, u0::Vector{<:Real},
                        no_AB_period::Real, AB_period::Real,
                        params::SimulationParams)

    all_ts = Float64[]
    all_ys = Vector{Vector{Float64}}()

    current_t = float(tspan[1])
    current_y = Float64.(u0)

    while current_t < tspan[2]
        end_t = min(current_t + no_AB_period, tspan[2])
        current_y[9] = 0.0
        prob = ODEProblem(vector_field!, current_y, (current_t, end_t), params)
        sol = solve(prob, Tsit5(); adaptive=false, dt=0.2)
        append!(all_ts, sol.t)
        append!(all_ys, sol.u)
        current_y .= sol.u[end]
        current_t = sol.t[end]
        if current_t >= tspan[2]
            break
        end

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

function single_simulation(u0::Vector{<:Real},
                           tspan::Tuple{Real,Real},
                           params::SimulationParams,
                           no_AB::Real, AB::Real;
                           output_dir::AbstractString="output/ode",
                           save_csv::Bool=true,
                           save_plots::Bool=true)

    mkpath(output_dir)

    results = run_simulation(tspan, u0, no_AB, AB, params)

    if save_csv
        CSV.write(joinpath(output_dir, "simulation.csv"), results)
    end

    if save_plots
        p_full = plot(results.time,
            [results.nS0 results.nSR results.nSR2 results.nSS results.nR0 results.nRS results.nRR results.nRR2],
            xlabel="Time",
            ylabel="Population density",
            title="ODE simulation",
            label=["S0" "SR" "SR2" "SS" "R0" "RS" "RR" "RR2"],
            linewidth=2,
            legend=:outerright)

        p_log = plot(results.time,
            [results.nS0 results.nSR results.nSR2 results.nSS results.nR0 results.nRS results.nRR results.nRR2],
            xlabel="Time",
            ylabel="Population density (log10)",
            title="ODE simulation (log scale)",
            label=["S0" "SR" "SR2" "SS" "R0" "RS" "RR" "RR2"],
            yscale=:log10,
            ylims=(1e-12, 1.0),
            linewidth=2,
            legend=:outerright)

        savefig(p_full, joinpath(output_dir, "full.png"))
        savefig(p_log,  joinpath(output_dir, "full_log.png"))
    end

    return results
end

end
