using DifferentialEquations, DataFrames, Plots, Statistics

function vector_field!(du, u, p, t)

    nS0, nSR1, nSR2, nSS, nR0, nRS, nRR1, nRR2, A = u
    λ, α, β1, β2, s1, s2, x, cR = p
    N = sum(u[1:8])

    λP1   = (1 - x) * λ              # nSS: plasmid carriage cost
    λRP1  = (1 - cR) * (1 - x) * λ  # nSR1: resistance + plasmid carriage cost
    λP2   = (1 - x) * λ              # plasmid 2 carriage cost
    λRP2  = (1 - cR) * (1 - x) * λ  # nSR2: resistance + plasmid carriage cost
    λR    = (1 - cR) * λ             # nR0: chromosomal resistance cost only
    λRRP1 = (1 - cR) * (1 - x) * λ  ## No double cost
    λRRP2 = (1 - cR) * (1 - x) * λ  ## No double cost

    donors1 = nSS + nSR1 + nRS + nRR1
    donors2 = nSR2 + nRR2

    ## nS0
    du[1] = nS0 * (λ - β1*donors1 - β2*donors2 - α*N - A) +
            s1*(λP1*nSS + λRP1*nSR1) + s2*λRP2*nSR2

    ## nSR1
    du[2] = nSR1 * (λRP1*(1 - s1) - α*N) + β1*(nSR1 + nRR1)*nS0

    ## nSR2
    du[3] = nSR2 * (λRP2*(1 - s2) - α*N) + β2*(nSR2 + nRR2)*nS0

    ## nSS
    du[4] = nSS * (λP1*(1 - s1) - α*N - A) + β1*(nSS + nRS)*nS0

    ## nR0
    du[5] = nR0 * (λR - β1*donors1 - β2*donors2 - α*N) +
            s1*(λRRP1*nRS + λRRP1*nRR1) + s2*λRRP2*nRR2

    ## nRS
    du[6] = nRS * (λRRP1*(1 - s1) - α*N) + β1*(nSS + nRS)*nR0

    ## nRR1
    du[7] = nRR1 * (λRRP1*(1 - s1) - α*N) + β1*(nSR1 + nRR1)*nR0

    ## nRR2
    du[8] = nRR2 * (λRRP2*(1 - s2) - α*N) + β2*(nSR2 + nRR2)*nR0

end

function run_simulation(tspan, u0, no_AB_period, AB_period, args)

    all_ts = Vector{Float64}()
    all_ys = Vector{Vector{Float64}}()

    current_t = tspan[1]
    current_y = copy(u0)

    while current_t < tspan[2]

        ## No AB period
        end_t = min(current_t + no_AB_period, tspan[2])
        current_y[end] = 0.0

        prob = ODEProblem(vector_field!, current_y, (current_t, end_t), args)
        sol  = solve(prob, Tsit5(); adaptive = false, dt = 0.2)

        append!(all_ts, sol.t)
        append!(all_ys, sol.u)

        current_y = sol.u[end]
        current_t = sol.t[end]
        current_t >= tspan[2] && break

        ## AB period
        end_t = min(current_t + AB_period, tspan[2])
        current_y[end] = 1.0

        prob = ODEProblem(vector_field!, current_y, (current_t, end_t), args)
        sol  = solve(prob, Tsit5(); adaptive = false, dt = 0.2)

        append!(all_ts, sol.t)
        append!(all_ys, sol.u)

        current_y = sol.u[end]
        current_t = sol.t[end]

    end

    col_names = [:time, :nS0, :nSR1, :nSR2, :nSS, :nR0, :nRS, :nRR1, :nRR2, :A]
    data      = hcat(all_ts, reduce(hcat, all_ys)')
    return DataFrame(data, col_names)

end