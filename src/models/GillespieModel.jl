## gillespie_standalone.jl
using Random, Distributions, CSV, DataFrames

struct ModelParameters
    λ::Float64
    λSS::Float64
    λR::Float64
    λRR::Float64
    λRR2::Float64

    SS_s::Float64
    SS_β::Float64
    s1::Float64
    β1::Float64
    s2::Float64
    β2::Float64

    α::Float64
    dose::Float64
    K::Float64

    AB_period::Float64
    no_AB_period::Float64
end

@inline function get_A(absolute_time::Float64, dose::Float64, AB_period::Float64, no_AB_period::Float64)
    AB_period == 0.0 && return 0.0
    cycle_len = AB_period + no_AB_period
    pos = mod(absolute_time, cycle_len)
    return pos >= no_AB_period ? dose : 0.0
end

@inline function time_to_next_ab_switch(absolute_time::Float64, AB_period::Float64, no_AB_period::Float64)
    AB_period == 0.0 && return Inf
    cycle_len = AB_period + no_AB_period
    pos = mod(absolute_time, cycle_len)
    return pos < no_AB_period ? (no_AB_period - pos) : (cycle_len - pos)
end

# x = [S0, SR, SR2, SS, R0, RS, RR, RR2]
function gillespie_patch!(x::Vector{Int32}, p::ModelParameters, tf::Float64,
                         model_time::Float64, rng::AbstractRNG)

    λ, λSS, λR, λRR, λRR2 = p.λ, p.λSS, p.λR, p.λRR, p.λRR2
    SS_s, SS_β, s1, β1, s2, β2 = p.SS_s, p.SS_β, p.s1, p.β1, p.s2, p.β2
    α, K, dose = p.α, p.K, p.dose
    AB_period, no_AB_period = p.AB_period, p.no_AB_period

    nS0, nSR, nSR2, nSS_, nR0, nRS, nRR_, nRR2_ = x

    dt = 0.0
    invK = 1.0 / K
    βSS_over_K = SS_β * invK
    β1_over_K  = β1  * invK
    β2_over_K  = β2  * invK

    while dt < tf
        N = nS0 + nSR + nSR2 + nSS_ + nR0 + nRS + nRR_ + nRR2_
        abs_time = model_time + dt
        A = get_A(abs_time, dose, AB_period, no_AB_period)

        αN_over_K = α * N * invK

        nSR_plus_nRR   = nSR + nRR_
        nSR2_plus_nRR2 = nSR2 + nRR2_
        nSS_plus_nRS   = nSS_ + nRS

        total_rate = (
            λ*nS0 + (αN_over_K + A)*nS0 +
            β1_over_K*nS0*nSR_plus_nRR + β2_over_K*nS0*nSR2_plus_nRR2 + βSS_over_K*nS0*nSS_plus_nRS +
            λRR*nSR + αN_over_K*nSR + s1*nSR*λRR +
            λRR2*nSR2 + αN_over_K*nSR2 + s2*nSR2*λRR2 +
            λSS*nSS_ + (αN_over_K + A)*nSS_ + SS_s*nSS_*λSS +
            λR*nR0 + αN_over_K*nR0 +
            β1_over_K*nR0*nSR_plus_nRR + β2_over_K*nR0*nSR2_plus_nRR2 + βSS_over_K*nR0*nSS_plus_nRS +
            λRR*nRS + αN_over_K*nRS + SS_s*nRS*λRR +
            λRR*nRR_ + αN_over_K*nRR_ + s1*nRR_*λRR +
            λRR2*nRR2_ + αN_over_K*nRR2_ + s2*nRR2_*λRR2
        )

        total_rate <= 0.0 && break

        τ = -log(rand(rng)) / total_rate
        t_switch = time_to_next_ab_switch(abs_time, AB_period, no_AB_period)
        if τ > t_switch && (dt + t_switch) < tf
            dt += t_switch
            continue
        end

        dt += τ
        dt > tf && break

        ru = rand(rng) * total_rate
        cumulative = 0.0

        cumulative += λ*nS0
        if ru <= cumulative
            nS0 += 1; continue
        end
        cumulative += (αN_over_K + A)*nS0
        if ru <= cumulative
            nS0 = max(nS0 - 1, 0); continue
        end
        cumulative += β1_over_K*nS0*nSR_plus_nRR
        if ru <= cumulative
            nS0 = max(nS0 - 1, 0); nSR += 1; continue
        end
        cumulative += β2_over_K*nS0*nSR2_plus_nRR2
        if ru <= cumulative
            nS0 = max(nS0 - 1, 0); nSR2 += 1; continue
        end
        cumulative += βSS_over_K*nS0*nSS_plus_nRS
        if ru <= cumulative
            nS0 = max(nS0 - 1, 0); nSS_ += 1; continue
        end

        cumulative += λRR*nSR
        if ru <= cumulative
            nSR += 1; continue
        end
        cumulative += αN_over_K*nSR
        if ru <= cumulative
            nSR = max(nSR - 1, 0); continue
        end
        cumulative += s1*nSR*λRR
        if ru <= cumulative
            nSR = max(nSR - 1, 0); nS0 += 1; continue
        end

        cumulative += λRR2*nSR2
        if ru <= cumulative
            nSR2 += 1; continue
        end
        cumulative += αN_over_K*nSR2
        if ru <= cumulative
            nSR2 = max(nSR2 - 1, 0); continue
        end
        cumulative += s2*nSR2*λRR2
        if ru <= cumulative
            nSR2 = max(nSR2 - 1, 0); nS0 += 1; continue
        end

        cumulative += λSS*nSS_
        if ru <= cumulative
            nSS_ += 1; continue
        end
        cumulative += (αN_over_K + A)*nSS_
        if ru <= cumulative
            nSS_ = max(nSS_ - 1, 0); continue
        end
        cumulative += SS_s*nSS_*λSS
        if ru <= cumulative
            nSS_ = max(nSS_ - 1, 0); nS0 += 1; continue
        end

        cumulative += λR*nR0
        if ru <= cumulative
            nR0 += 1; continue
        end
        cumulative += αN_over_K*nR0
        if ru <= cumulative
            nR0 = max(nR0 - 1, 0); continue
        end
        cumulative += β1_over_K*nR0*nSR_plus_nRR
        if ru <= cumulative
            nR0 = max(nR0 - 1, 0); nRR_ += 1; continue
        end
        cumulative += β2_over_K*nR0*nSR2_plus_nRR2
        if ru <= cumulative
            nR0 = max(nR0 - 1, 0); nRR2_ += 1; continue
        end
        cumulative += βSS_over_K*nR0*nSS_plus_nRS
        if ru <= cumulative
            nR0 = max(nR0 - 1, 0); nRS += 1; continue
        end

        cumulative += λRR*nRS
        if ru <= cumulative
            nRS += 1; continue
        end
        cumulative += αN_over_K*nRS
        if ru <= cumulative
            nRS = max(nRS - 1, 0); continue
        end
        cumulative += SS_s*nRS*λRR
        if ru <= cumulative
            nRS = max(nRS - 1, 0); nR0 += 1; continue
        end

        cumulative += λRR*nRR_
        if ru <= cumulative
            nRR_ += 1; continue
        end
        cumulative += αN_over_K*nRR_
        if ru <= cumulative
            nRR_ = max(nRR_ - 1, 0); continue
        end
        cumulative += s1*nRR_*λRR
        if ru <= cumulative
            nRR_ = max(nRR_ - 1, 0); nR0 += 1; continue
        end

        cumulative += λRR2*nRR2_
        if ru <= cumulative
            nRR2_ += 1; continue
        end
        cumulative += αN_over_K*nRR2_
        if ru <= cumulative
            nRR2_ = max(nRR2_ - 1, 0); continue
        end
        nRR2_ = max(nRR2_ - 1, 0); nR0 += 1
    end

    x[1]=nS0; x[2]=nSR; x[3]=nSR2; x[4]=nSS_
    x[5]=nR0; x[6]=nRS; x[7]=nRR_; x[8]=nRR2_
    return x
end

function run_patch(x0::Vector{Int32}, p::ModelParameters;
                   dt::Float64=1.0, t_end::Float64=500.0, seed::Int=3,
                   output_dir::AbstractString="output_gillespie", save_csv::Bool=true)
    isdir(output_dir) || mkpath(output_dir)

    rng = MersenneTwister(seed)
    x = copy(x0)

    rows = Vector{NTuple{9,Float64}}()
    t = 0.0
    while t < t_end
        push!(rows, (t, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]))
        gillespie_patch!(x, p, dt, t, rng)
        t += dt
    end

    df = DataFrame(rows, [:time,:nS0,:nSR,:nSR2,:nSS,:nR0,:nRS,:nRR,:nRR2])
    save_csv && CSV.write(joinpath(output_dir, "gillespie_patch.csv"), df)
    return df
end

# =========================
# MAIN
# =========================
println("\n" * "="^70)
println("GILLESPIE (SINGLE PATCH): RUN")
println("="^70)

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
    10000.0,# K
    40.0,       # AB_period
    350.0       # no_AB_period
)

x0 = Int32[0, 1000, 0, 0, 1000, 0, 0, 0]  # [S0, SR, SR2, SS, R0, RS, RR, RR2]
@time df = run_patch(x0, params; dt=1.0, t_end=500.0, seed=3, output_dir="output_gillespie", save_csv=true)
println("Done. Rows: ", nrow(df))
