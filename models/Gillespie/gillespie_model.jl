using Random, Statistics, Distributions, DataFrames

## Parameters
struct GillespieParameters
    λ::Float64      # base growth rate (nS0 cells)
    λP::Float64     # growth rate with sensitive plasmid (nSS, nRS)
    λR::Float64     # growth rate with chrom. resistance, no plasmid (nR0)
    λRP::Float64    # growth rate with resistance + plasmid (nSR1, nSR2)
    λRRP::Float64   # growth rate with chrom. resistance + plasmid (nRS, nRR1, nRR2); no double cost

    s_ss::Float64   # segregation rate, sensitive plasmid
    β_ss::Float64   # conjugation rate, sensitive plasmid
    s1::Float64     # segregation rate, resistant plasmid 1
    β1::Float64     # conjugation rate, resistant plasmid 1
    s2::Float64     # segregation rate, resistant plasmid 2
    β2::Float64     # conjugation rate, resistant plasmid 2

    α::Float64      # density-dependent death rate coefficient
    dose::Float64   # antibiotic dose (kills nS0 and nSS)
    K::Float64      # carrying capacity (scales death and conjugation)

    AB_period::Float64
    no_AB_period::Float64
end

@inline function get_A(t::Float64, phase::Int8, dose::Float64,
                       AB_period::Float64, no_AB_period::Float64)
    AB_period == 0.0 && return 0.0
    cycle = AB_period + no_AB_period
    pos   = mod(t, cycle)
    phase == 0 ? (pos >= no_AB_period ? dose : 0.0) :
                 (pos < AB_period     ? dose : 0.0)
end

@inline function time_to_next_switch(t::Float64, phase::Int8,
                                     AB_period::Float64, no_AB_period::Float64)
    AB_period == 0.0 && return Inf
    cycle = AB_period + no_AB_period
    pos   = mod(t, cycle)
    phase == 0 ? (pos < no_AB_period ? no_AB_period - pos : cycle - pos) :
                 (pos < AB_period    ? AB_period - pos    : cycle - pos)
end

function gillespie!(x::Vector{Int32}, params::GillespieParameters, tf::Float64,
                    t0::Float64, phase::Int8, rng::AbstractRNG)

    λ, λP, λR, λRP, λRRP = params.λ, params.λP, params.λR, params.λRP, params.λRRP
    s_ss, β_ss = params.s_ss, params.β_ss
    s1, β1     = params.s1, params.β1
    s2, β2     = params.s2, params.β2
    α, K, dose = params.α, params.K, params.dose
    AB_period  = params.AB_period
    no_AB_period = params.no_AB_period

    nS0, nSR1, nSR2, nSS, nR0, nRS, nRR1, nRR2 = x

    invK      = 1.0 / K
    β_ss_K    = β_ss * invK
    β1_K      = β1  * invK
    β2_K      = β2  * invK

    dt = 0.0

    while dt < tf

        N   = nS0 + nSR1 + nSR2 + nSS + nR0 + nRS + nRR1 + nRR2
        A   = get_A(t0 + dt, phase, dose, AB_period, no_AB_period)
        αNK = α * N * invK

        donors_r1 = nSR1 + nRR1
        donors_r2 = nSR2 + nRR2
        donors_ss = nSS  + nRS

        total_rate = (
            ## nS0
            λ*nS0 + (αNK + A)*nS0 +
            β1_K*nS0*donors_r1 + β2_K*nS0*donors_r2 + β_ss_K*nS0*donors_ss +
            ## nSR1
            λRP*nSR1 + αNK*nSR1 + s1*nSR1*λRP +
            ## nSR2
            λRP*nSR2 + αNK*nSR2 + s2*nSR2*λRP +
            ## nSS
            λP*nSS + (αNK + A)*nSS + s_ss*nSS*λP +
            ## nR0
            λR*nR0 + αNK*nR0 +
            β1_K*nR0*donors_r1 + β2_K*nR0*donors_r2 + β_ss_K*nR0*donors_ss +
            ## nRS
            λRRP*nRS + αNK*nRS + s_ss*nRS*λRRP +
            ## nRR1
            λRRP*nRR1 + αNK*nRR1 + s1*nRR1*λRRP +
            ## nRR2
            λRRP*nRR2 + αNK*nRR2 + s2*nRR2*λRRP
        )

        total_rate <= 0.0 && break

        τ = -log(rand(rng)) / total_rate

        ## If an AB switch happens before the next event, skip ahead to it
        t_switch = time_to_next_switch(t0 + dt, phase, AB_period, no_AB_period)
        if τ > t_switch && (dt + t_switch) < tf
            dt += t_switch
            continue
        end

        dt += τ
        dt > tf && break

        ## Select event
        ru  = rand(rng) * total_rate
        cumulative = 0.0

        ## nS0 births
        cumulative += λ*nS0;              ru <= cumulative && (nS0 += 1;  continue)
        ## nS0 deaths (density-dependent + AB)
        cumulative += (αNK + A)*nS0;      ru <= cumulative && (nS0 = max(0, nS0-1); continue)
        ## nS0 acquires resistant plasmid 1 → becomes nSR1
        cumulative += β1_K*nS0*donors_r1; ru <= cumulative && (nS0 = max(0, nS0-1); nSR1 += 1; continue)
        ## nS0 acquires resistant plasmid 2 → becomes nSR2
        cumulative += β2_K*nS0*donors_r2; ru <= cumulative && (nS0 = max(0, nS0-1); nSR2 += 1; continue)
        ## nS0 acquires sensitive plasmid → becomes nSS
        cumulative += β_ss_K*nS0*donors_ss; ru <= cumulative && (nS0 = max(0, nS0-1); nSS += 1; continue)

        ## nSR1 births
        cumulative += λRP*nSR1;          ru <= cumulative && (nSR1 += 1; continue)
        ## nSR1 deaths
        cumulative += αNK*nSR1;           ru <= cumulative && (nSR1 = max(0, nSR1-1); continue)
        ## nSR1 segregates plasmid → becomes nS0
        cumulative += s1*nSR1*λRP;       ru <= cumulative && (nSR1 = max(0, nSR1-1); nS0 += 1; continue)

        ## nSR2 births
        cumulative += λRP*nSR2;          ru <= cumulative && (nSR2 += 1; continue)
        ## nSR2 deaths
        cumulative += αNK*nSR2;           ru <= cumulative && (nSR2 = max(0, nSR2-1); continue)
        ## nSR2 segregates plasmid → becomes nS0
        cumulative += s2*nSR2*λRP;       ru <= cumulative && (nSR2 = max(0, nSR2-1); nS0 += 1; continue)

        ## nSS births
        cumulative += λP*nSS;            ru <= cumulative && (nSS += 1; continue)
        ## nSS deaths (density-dependent + AB)
        cumulative += (αNK + A)*nSS;      ru <= cumulative && (nSS = max(0, nSS-1); continue)
        ## nSS segregates plasmid → becomes nS0
        cumulative += s_ss*nSS*λP;       ru <= cumulative && (nSS = max(0, nSS-1); nS0 += 1; continue)

        ## nR0 births
        cumulative += λR*nR0;            ru <= cumulative && (nR0 += 1; continue)
        ## nR0 deaths
        cumulative += αNK*nR0;            ru <= cumulative && (nR0 = max(0, nR0-1); continue)
        ## nR0 acquires resistant plasmid 1 → becomes nRR1
        cumulative += β1_K*nR0*donors_r1; ru <= cumulative && (nR0 = max(0, nR0-1); nRR1 += 1; continue)
        ## nR0 acquires resistant plasmid 2 → becomes nRR2
        cumulative += β2_K*nR0*donors_r2; ru <= cumulative && (nR0 = max(0, nR0-1); nRR2 += 1; continue)
        ## nR0 acquires sensitive plasmid → becomes nRS
        cumulative += β_ss_K*nR0*donors_ss; ru <= cumulative && (nR0 = max(0, nR0-1); nRS += 1; continue)

        ## nRS births
        cumulative += λRRP*nRS;           ru <= cumulative && (nRS += 1; continue)
        ## nRS deaths
        cumulative += αNK*nRS;            ru <= cumulative && (nRS = max(0, nRS-1); continue)
        ## nRS segregates plasmid → becomes nR0
        cumulative += s_ss*nRS*λRRP;      ru <= cumulative && (nRS = max(0, nRS-1); nR0 += 1; continue)

        ## nRR1 births
        cumulative += λRRP*nRR1;          ru <= cumulative && (nRR1 += 1; continue)
        ## nRR1 deaths
        cumulative += αNK*nRR1;           ru <= cumulative && (nRR1 = max(0, nRR1-1); continue)
        ## nRR1 segregates plasmid → becomes nR0
        cumulative += s1*nRR1*λRRP;       ru <= cumulative && (nRR1 = max(0, nRR1-1); nR0 += 1; continue)

        ## nRR2 births
        cumulative += λRRP*nRR2;          ru <= cumulative && (nRR2 += 1; continue)
        ## nRR2 deaths
        cumulative += αNK*nRR2;           ru <= cumulative && (nRR2 = max(0, nRR2-1); continue)
        ## nRR2 segregates plasmid → becomes nR0 (final event, always selected if reached)
        nRR2 = max(0, nRR2-1); nR0 += 1

    end

    x[1]=nS0; x[2]=nSR1; x[3]=nSR2; x[4]=nSS
    x[5]=nR0; x[6]=nRS;  x[7]=nRR1; x[8]=nRR2

    return x
end

## Run

function run_gillespie(x0::Vector{Int32}, params::GillespieParameters, T::Float64;
                       record_interval::Float64 = 10.0,
                       phase::Int8 = Int8(0),
                       seed::Int   = 1234)

    rng    = MersenneTwister(seed)
    x      = copy(x0)
    times  = Float64[0.0]
    states = [copy(x)]

    t = 0.0
    while t < T
        step = min(record_interval, T - t)
        gillespie!(x, params, step, t, phase, rng)
        t += step
        push!(times, t)
        push!(states, copy(x))
        all(==(0), x) && (println("Extinction at t=$t"); break)
    end

    col_names = [:time, :nS0, :nSR1, :nSR2, :nSS, :nR0, :nRS, :nRR1, :nRR2]
    data = hcat(times, reduce(hcat, states)')
    return DataFrame(data, col_names)

end
