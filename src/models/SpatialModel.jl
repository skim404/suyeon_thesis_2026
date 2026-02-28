## spatial_standalone.jl
using CairoMakie, Colors, Random, Statistics, Distributions, Base.Threads

mutable struct GridSpace
    rows::Int
    cols::Int
    nS0::Matrix{Int32}
    nSR::Matrix{Int32}
    nSR2::Matrix{Int32}
    nSS::Matrix{Int32}
    nR0::Matrix{Int32}
    nRS::Matrix{Int32}
    nRR::Matrix{Int32}
    nRR2::Matrix{Int32}
    total_pop::Matrix{Int32}
    type_matrix::Matrix{UInt8}
    color_matrix::Matrix{RGBA{Float32}}
    normalized_intensity::Matrix{Float32}

    function GridSpace(rows::Int, cols::Int)
        new(
            rows, cols,
            zeros(Int32, rows, cols),
            zeros(Int32, rows, cols),
            zeros(Int32, rows, cols),
            zeros(Int32, rows, cols),
            zeros(Int32, rows, cols),
            zeros(Int32, rows, cols),
            zeros(Int32, rows, cols),
            zeros(Int32, rows, cols),
            zeros(Int32, rows, cols),
            zeros(UInt8, rows, cols),
            fill(RGBA{Float32}(1,1,1,0), rows, cols),
            zeros(Float32, rows, cols)
        )
    end
end

mutable struct TimeSeriesData
    times::Vector{Float64}
    nS0::Vector{Float64}
    nSR::Vector{Float64}
    nSR2::Vector{Float64}
    nSS::Vector{Float64}
    nR0::Vector{Float64}
    nRS::Vector{Float64}
    nRR::Vector{Float64}
    nRR2::Vector{Float64}
    function TimeSeriesData()
        new(Float64[], Float64[], Float64[], Float64[], Float64[],
            Float64[], Float64[], Float64[], Float64[])
    end
end

function update_time_series!(ts::TimeSeriesData, grid::GridSpace, time::Float64)
    total_patches = grid.rows * grid.cols
    push!(ts.times, time)
    push!(ts.nS0,  sum(grid.nS0)  / total_patches)
    push!(ts.nSR,  sum(grid.nSR)  / total_patches)
    push!(ts.nSR2, sum(grid.nSR2) / total_patches)
    push!(ts.nSS,  sum(grid.nSS)  / total_patches)
    push!(ts.nR0,  sum(grid.nR0)  / total_patches)
    push!(ts.nRS,  sum(grid.nRS)  / total_patches)
    push!(ts.nRR,  sum(grid.nRR)  / total_patches)
    push!(ts.nRR2, sum(grid.nRR2) / total_patches)
end

function update_totals!(grid::GridSpace)
    @inbounds for j in 1:grid.cols, i in 1:grid.rows
        grid.total_pop[i,j] =
            grid.nS0[i,j] + grid.nSR[i,j] + grid.nSR2[i,j] + grid.nSS[i,j] +
            grid.nR0[i,j] + grid.nRS[i,j] + grid.nRR[i,j] + grid.nRR2[i,j]
    end
end

function initialize_grid!(grid::GridSpace;
                         population_fractions::Vector{Float64},
                         population_types::Vector{Vector{Int}},
                         seed::Int=1)

    rng = MersenneTwister(seed)
    n_patches = grid.rows * grid.cols
    assignments = zeros(Int, n_patches)
    idx = 1

    for (type_idx, frac) in enumerate(population_fractions)
        n_this = round(Int, frac * n_patches)
        for _ in 1:n_this
            idx <= n_patches || break
            assignments[idx] = type_idx
            idx += 1
        end
    end

    Random.shuffle!(rng, assignments)

    @inbounds for lin in 1:n_patches
        i = div(lin - 1, grid.cols) + 1
        j = mod(lin - 1, grid.cols) + 1
        a = assignments[lin]
        if a > 0
            pop = population_types[a]
            grid.nS0[i,j]  = pop[1]
            grid.nSR[i,j]  = pop[2]
            grid.nSR2[i,j] = pop[3]
            grid.nSS[i,j]  = pop[4]
            grid.nR0[i,j]  = pop[5]
            grid.nRS[i,j]  = pop[6]
            grid.nRR[i,j]  = pop[7]
            grid.nRR2[i,j] = pop[8]
        end
    end

    update_totals!(grid)
end

function get_dominant_types!(grid::GridSpace)
    @inbounds for j in 1:grid.cols, i in 1:grid.rows
        if grid.total_pop[i,j] == 0
            grid.type_matrix[i,j] = 0
            continue
        end

        p1 = grid.nSR[i,j]
        p2 = grid.nSR2[i,j]
        chr = grid.nR0[i,j] + grid.nRS[i,j] + grid.nRR[i,j] + grid.nRR2[i,j]
        sens = grid.nS0[i,j] + grid.nSS[i,j]

        m = max(p1, p2, chr, sens)
        grid.type_matrix[i,j] = if p1 == m && p1 > 0
            1
        elseif p2 == m && p2 > 0
            2
        elseif chr == m && chr > 0
            3
        else
            4
        end
    end
end

function create_color_matrix!(grid::GridSpace)
    max_pop = maximum(grid.total_pop)
    grid.normalized_intensity .= max_pop > 0 ? (grid.total_pop ./ max_pop) : 0.0

    @inbounds for j in 1:grid.cols, i in 1:grid.rows
        t = grid.type_matrix[i,j]
        op = grid.normalized_intensity[i,j]
        if t == 0
            grid.color_matrix[i,j] = RGBA(1,1,1,0)
        elseif t == 1
            grid.color_matrix[i,j] = RGBA(0,0,1,op)
        elseif t == 2
            grid.color_matrix[i,j] = RGBA(0,1,0,op)
        elseif t == 3
            gcomp = 1.0 - 0.5 * op
            grid.color_matrix[i,j] = RGBA(1,gcomp,0,op)
        else
            grid.color_matrix[i,j] = RGBA(0.3,0.3,0.3,op)
        end
    end
end

function setup_visualization(grid::GridSpace, ts::TimeSeriesData)
    get_dominant_types!(grid)
    create_color_matrix!(grid)

    fig = Figure(size=(1200, 500))

    ax1 = Axis(fig[1, 1], aspect=DataAspect(),
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false,
               xgridvisible=false, ygridvisible=false,
               leftspinevisible=false, rightspinevisible=false,
               topspinevisible=false, bottomspinevisible=false,
               title="Spatial distribution")

    implot = image!(ax1, grid.color_matrix, interpolate=false)

    ax2 = Axis(fig[1, 2], xlabel="Time", ylabel="Avg density per patch",
               title="Population dynamics")

    lines!(ax2, ts.times, ts.nS0, label="S0")
    lines!(ax2, ts.times, ts.nSR, label="SR")
    lines!(ax2, ts.times, ts.nSR2, label="SR2")
    lines!(ax2, ts.times, ts.nSS, label="SS")
    lines!(ax2, ts.times, ts.nR0, label="R0")
    lines!(ax2, ts.times, ts.nRS, label="RS")
    lines!(ax2, ts.times, ts.nRR, label="RR")
    lines!(ax2, ts.times, ts.nRR2, label="RR2")

    axislegend(ax2, position=:lt, labelsize=10)
    return fig, implot
end

struct ModelParameters
    λ::Float64
    λSS::Float64
    λR::Float64
    λRR::Float64
    λRR2::Float64

    SS_s::Float64; SS_β::Float64
    s1::Float64; β1::Float64
    s2::Float64; β2::Float64

    α::Float64
    dose::Float64
    K::Float64

    migration_rate::Float64
    extinction_rate::Float64
    dt::Float64

    AB_period::Float64
    no_AB_period::Float64
    out_of_phase::Bool
    phase_fraction::Float64
end

mutable struct GridDynamics
    grid::GridSpace
    params::ModelParameters
    model_time::Float64
    ab_phase::Matrix{Int8}
    rngs::Vector{MersenneTwister}
    mig_rng::MersenneTwister
    thread_states::Vector{Vector{Int32}}
    time_series::TimeSeriesData
end

function init_ab_phase(rows::Int, cols::Int, params::ModelParameters; seed::Int=1)
    A = fill(Int8(0), rows, cols)
    if params.out_of_phase && params.phase_fraction > 0
        rng = MersenneTwister(seed)
        @inbounds for j in 1:cols, i in 1:rows
            A[i,j] = rand(rng) < params.phase_fraction ? Int8(1) : Int8(0)
        end
    end
    return A
end

function GridDynamics(grid::GridSpace, params::ModelParameters; seed::Int=1234)
    rows, cols = grid.rows, grid.cols
    seeds = rand(MersenneTwister(seed), UInt, max(1, nthreads()))
    rngs = [MersenneTwister(s) for s in seeds]
    ab_phase = init_ab_phase(rows, cols, params; seed=seed+1)
    mig_rng  = MersenneTwister(UInt(seed+2))
    states = [Vector{Int32}(undef, 8) for _ in 1:nthreads()]
    ts = TimeSeriesData()
    return GridDynamics(grid, params, 0.0, ab_phase, rngs, mig_rng, states, ts)
end

@inline function get_A_patch(absolute_time::Float64, phase::Int8,
                             dose::Float64, AB_period::Float64, no_AB_period::Float64)
    AB_period == 0.0 && return 0.0
    cycle_len = AB_period + no_AB_period
    pos = mod(absolute_time, cycle_len)
    if phase == 0
        return pos >= no_AB_period ? dose : 0.0
    else
        return pos < AB_period ? dose : 0.0
    end
end

@inline function time_to_next_ab_switch_patch(absolute_time::Float64, phase::Int8,
                                              AB_period::Float64, no_AB_period::Float64)
    AB_period == 0.0 && return Inf
    cycle_len = AB_period + no_AB_period
    pos = mod(absolute_time, cycle_len)
    if phase == 0
        return pos < no_AB_period ? (no_AB_period - pos) : (cycle_len - pos)
    else
        return pos < AB_period ? (AB_period - pos) : (cycle_len - pos)
    end
end

# buf = [S0, SR, SR2, SS, R0, RS, RR, RR2]
function gillespie!(x::Vector{Int32}, p::ModelParameters, tf::Float64,
                    model_time::Float64, ab_phase::Int8, rng::AbstractRNG)

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
        A = get_A_patch(abs_time, ab_phase, dose, AB_period, no_AB_period)

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
        t_switch = time_to_next_ab_switch_patch(abs_time, ab_phase, AB_period, no_AB_period)
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

@inline function handle_extinctions!(sim::GridDynamics)
    p_ext = sim.params.extinction_rate * sim.params.dt
    p_ext <= 0.0 && return
    g = sim.grid
    rng = sim.mig_rng
    @inbounds for j in 1:g.cols, i in 1:g.rows
        if g.total_pop[i,j] > 0 && rand(rng) < p_ext
            g.nS0[i,j]=0; g.nSR[i,j]=0; g.nSR2[i,j]=0; g.nSS[i,j]=0
            g.nR0[i,j]=0; g.nRS[i,j]=0; g.nRR[i,j]=0; g.nRR2[i,j]=0
            g.total_pop[i,j]=0
        end
    end
end

@inline lin_idx(i, j, cols) = (i - 1) * cols + j
@inline cart_idx(idx, cols) = (div(idx - 1, cols) + 1, mod(idx - 1, cols) + 1)

function handle_migration!(sim::GridDynamics)
    g = sim.grid; p = sim.params
    r = g.rows; c = g.cols
    p.migration_rate <= 0.0 && return

    rng = sim.mig_rng
    events = Vector{NTuple{10,Int32}}()
    sizehint!(events, r*c*2)

    @inbounds for j in 1:c, i in 1:r
        if g.total_pop[i,j] == 0
            continue
        end

        nbrs = (
            (i, mod1(j+1, c)),
            (i, mod1(j-1, c)),
            (mod1(i+1, r), j),
            (mod1(i-1, r), j)
        )

        nS0=g.nS0[i,j]; nSR=g.nSR[i,j]; nSR2=g.nSR2[i,j]; nSS=g.nSS[i,j]
        nR0=g.nR0[i,j]; nRS=g.nRS[i,j]; nRR=g.nRR[i,j]; nRR2=g.nRR2[i,j]

        total_prob = p.migration_rate * p.dt

        mS0  = nS0  > 0 ? min(rand(rng, Poisson(nS0  * total_prob)), nS0)   : 0
        mSR  = nSR  > 0 ? min(rand(rng, Poisson(nSR  * total_prob)), nSR)   : 0
        mSR2 = nSR2 > 0 ? min(rand(rng, Poisson(nSR2 * total_prob)), nSR2)  : 0
        mSS  = nSS  > 0 ? min(rand(rng, Poisson(nSS  * total_prob)), nSS)   : 0
        mR0  = nR0  > 0 ? min(rand(rng, Poisson(nR0  * total_prob)), nR0)   : 0
        mRS  = nRS  > 0 ? min(rand(rng, Poisson(nRS  * total_prob)), nRS)   : 0
        mRR  = nRR  > 0 ? min(rand(rng, Poisson(nRR  * total_prob)), nRR)   : 0
        mRR2 = nRR2 > 0 ? min(rand(rng, Poisson(nRR2 * total_prob)), nRR2)  : 0

        probs = [0.25,0.25,0.25,0.25]
        dS0  = mS0  > 0 ? rand(rng, Multinomial(mS0,  probs)) : [0,0,0,0]
        dSR  = mSR  > 0 ? rand(rng, Multinomial(mSR,  probs)) : [0,0,0,0]
        dSR2 = mSR2 > 0 ? rand(rng, Multinomial(mSR2, probs)) : [0,0,0,0]
        dSS  = mSS  > 0 ? rand(rng, Multinomial(mSS,  probs)) : [0,0,0,0]
        dR0  = mR0  > 0 ? rand(rng, Multinomial(mR0,  probs)) : [0,0,0,0]
        dRS  = mRS  > 0 ? rand(rng, Multinomial(mRS,  probs)) : [0,0,0,0]
        dRR  = mRR  > 0 ? rand(rng, Multinomial(mRR,  probs)) : [0,0,0,0]
        dRR2 = mRR2 > 0 ? rand(rng, Multinomial(mRR2, probs)) : [0,0,0,0]

        for (k,(ni,nj)) in enumerate(nbrs)
            m1=dS0[k]; m2=dSR[k]; m3=dSR2[k]; m4=dSS[k]; m5=dR0[k]; m6=dRS[k]; m7=dRR[k]; m8=dRR2[k]
            tot = m1+m2+m3+m4+m5+m6+m7+m8
            if tot > 0
                push!(events, (Int32(lin_idx(i,j,c)), Int32(lin_idx(ni,nj,c)),
                               Int32(m1),Int32(m2),Int32(m3),Int32(m4),
                               Int32(m5),Int32(m6),Int32(m7),Int32(m8)))
            end
        end
    end

    @inbounds for ev in events
        s_lin = Int(ev[1]); t_lin = Int(ev[2])
        si,sj = cart_idx(s_lin, c); ti,tj = cart_idx(t_lin, c)

        mS0  = min(Int(ev[3]),  g.nS0[si,sj]);  g.nS0[si,sj]  -= mS0;  g.nS0[ti,tj]  += mS0
        mSR  = min(Int(ev[4]),  g.nSR[si,sj]);  g.nSR[si,sj]  -= mSR;  g.nSR[ti,tj]  += mSR
        mSR2 = min(Int(ev[5]),  g.nSR2[si,sj]); g.nSR2[si,sj] -= mSR2; g.nSR2[ti,tj] += mSR2
        mSS  = min(Int(ev[6]),  g.nSS[si,sj]);  g.nSS[si,sj]  -= mSS;  g.nSS[ti,tj]  += mSS
        mR0  = min(Int(ev[7]),  g.nR0[si,sj]);  g.nR0[si,sj]  -= mR0;  g.nR0[ti,tj]  += mR0
        mRS  = min(Int(ev[8]),  g.nRS[si,sj]);  g.nRS[si,sj]  -= mRS;  g.nRS[ti,tj]  += mRS
        mRR  = min(Int(ev[9]),  g.nRR[si,sj]);  g.nRR[si,sj]  -= mRR;  g.nRR[ti,tj]  += mRR
        mRR2 = min(Int(ev[10]), g.nRR2[si,sj]); g.nRR2[si,sj] -= mRR2; g.nRR2[ti,tj] += mRR2
    end

    update_totals!(g)
end

function step!(sim::GridDynamics)
    g, p = sim.grid, sim.params
    rows, cols = g.rows, g.cols

    @threads :static for j in 1:cols
        rng = sim.rngs[threadid()]
        buf = sim.thread_states[threadid()]
        @inbounds for i in 1:rows
            if g.total_pop[i,j] == 0
                continue
            end

            buf[1]=g.nS0[i,j];  buf[2]=g.nSR[i,j];  buf[3]=g.nSR2[i,j]; buf[4]=g.nSS[i,j]
            buf[5]=g.nR0[i,j];  buf[6]=g.nRS[i,j];  buf[7]=g.nRR[i,j];  buf[8]=g.nRR2[i,j]

            gillespie!(buf, p, p.dt, sim.model_time, sim.ab_phase[i,j], rng)

            g.nS0[i,j]=max(0,buf[1]); g.nSR[i,j]=max(0,buf[2]); g.nSR2[i,j]=max(0,buf[3]); g.nSS[i,j]=max(0,buf[4])
            g.nR0[i,j]=max(0,buf[5]); g.nRS[i,j]=max(0,buf[6]); g.nRR[i,j]=max(0,buf[7]); g.nRR2[i,j]=max(0,buf[8])
        end
    end

    update_totals!(g)
    handle_extinctions!(sim)
    handle_migration!(sim)
    sim.model_time += p.dt
    return sim
end

is_extinct(grid::GridSpace) = all(==(0), grid.total_pop)

function run_simulation!(sim::GridDynamics; n_steps::Int=200, display_interval::Int=50,
                         output_dir::String="output_spatial_frames")

    isdir(output_dir) || mkpath(output_dir)
    ts = sim.time_series

    update_time_series!(ts, sim.grid, sim.model_time)
    fig, implot = setup_visualization(sim.grid, ts)

    function save_frame(stepnum::Int)
        get_dominant_types!(sim.grid)
        create_color_matrix!(sim.grid)
        implot[1] = sim.grid.color_matrix
        fname = joinpath(output_dir, "frame_$(lpad(stepnum, 4, '0'))_t_$(round(Int, sim.model_time)).png")
        save(fname, fig, px_per_unit=2)
        println("Saved: ", fname)
    end

    save_frame(0)

    for stepnum in 1:n_steps
        step!(sim)
        update_time_series!(ts, sim.grid, sim.model_time)

        if is_extinct(sim.grid)
            println("Extinction detected at step ", stepnum, " (t=", sim.model_time, ")")
            save_frame(stepnum)
            break
        end

        if display_interval > 0 && (stepnum % display_interval == 0)
            save_frame(stepnum)
        end
    end

    return sim
end

# =========================
# MAIN
# =========================
println("\n" * "="^70)
println("SPATIAL MODEL: RUN")
println("="^70)

master_seed = 3

grid = GridSpace(20, 20)
initialize_grid!(grid;
    population_fractions = [0.5, 0.25, 0.25],
    population_types = [
        [3000, 3000, 0,    0,    0,    0,    0, 0],   # S0 + SR
        [3000, 0,    3000, 0,    0,    0,    0, 0],   # S0 + SR2
        [0,    0,    0,    0, 3000, 3000,    0, 0]    # R0 + RS
    ],
    seed=master_seed
)

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
    0.001,      # migration_rate
    0.0,        # extinction_rate
    1.0,        # dt
    40.0,       # AB_period
    350.0,      # no_AB_period
    true,       # out_of_phase
    0.5         # phase_fraction
)

sim = GridDynamics(grid, params; seed=master_seed)
@time run_simulation!(sim; n_steps=200, display_interval=50, output_dir="output_spatial_frames")
println("Done.")
