using Random, Statistics, Base.Threads, Distributions, Colors

mutable struct GridSpace
    rows::Int
    cols::Int
    nS0::Matrix{Int32}
    nSR1::Matrix{Int32}
    nSR2::Matrix{Int32}
    nSS::Matrix{Int32}
    nR0::Matrix{Int32}
    nRS::Matrix{Int32}
    nRR1::Matrix{Int32}
    nRR2::Matrix{Int32}
    total_pop::Matrix{Int32}

    function GridSpace(rows::Int, cols::Int)
        new(rows, cols,
            zeros(Int32, rows, cols),  # nS0
            zeros(Int32, rows, cols),  # nSR1
            zeros(Int32, rows, cols),  # nSR2
            zeros(Int32, rows, cols),  # nSS
            zeros(Int32, rows, cols),  # nR0
            zeros(Int32, rows, cols),  # nRS
            zeros(Int32, rows, cols),  # nRR1
            zeros(Int32, rows, cols),  # nRR2
            zeros(Int32, rows, cols))  # total_pop
    end
end

mutable struct TimeSeriesData
    times::Vector{Float64}
    nS0::Vector{Float64}
    nSR1::Vector{Float64}
    nSR2::Vector{Float64}
    nSS::Vector{Float64}
    nR0::Vector{Float64}
    nRS::Vector{Float64}
    nRR1::Vector{Float64}
    nRR2::Vector{Float64}

    function TimeSeriesData()
        new((Vector{Float64}() for _ in 1:9)...)
    end
end

function update_time_series!(ts::TimeSeriesData, grid::GridSpace, time::Float64)
    n = grid.rows * grid.cols
    push!(ts.times, time)
    push!(ts.nS0,  sum(grid.nS0)  / n)
    push!(ts.nSR1, sum(grid.nSR1) / n)
    push!(ts.nSR2, sum(grid.nSR2) / n)
    push!(ts.nSS,  sum(grid.nSS)  / n)
    push!(ts.nR0,  sum(grid.nR0)  / n)
    push!(ts.nRS,  sum(grid.nRS)  / n)
    push!(ts.nRR1, sum(grid.nRR1) / n)
    push!(ts.nRR2, sum(grid.nRR2) / n)
end

function update_totals!(grid::GridSpace)
    @inbounds for j in 1:grid.cols, i in 1:grid.rows
        grid.total_pop[i,j] = grid.nS0[i,j] + grid.nSR1[i,j] + grid.nSR2[i,j] + grid.nSS[i,j] +
                               grid.nR0[i,j] + grid.nRS[i,j]  + grid.nRR1[i,j] + grid.nRR2[i,j]
    end
end

is_extinct(grid::GridSpace) = all(==(0), grid.total_pop)
function get_type_matrix(grid::GridSpace)
    out = zeros(Int, grid.rows, grid.cols)
    @inbounds for j in 1:grid.cols, i in 1:grid.rows
        if grid.total_pop[i,j] == 0
            out[i,j] = 0
            continue
        end
        ## All carriers of each resistance category
        p1_count    = grid.nSR1[i,j] + grid.nRR1[i,j]
        p2_count    = grid.nSR2[i,j] + grid.nRR2[i,j]
        chrom_count = grid.nR0[i,j]  + grid.nRS[i,j] + grid.nRR1[i,j] + grid.nRR2[i,j]

        if p1_count == 0 && p2_count == 0 && chrom_count == 0
            out[i,j] = 4  # sensitive only
        elseif p1_count >= p2_count && p1_count >= chrom_count
            out[i,j] = 1
        elseif p2_count >= p1_count && p2_count >= chrom_count
            out[i,j] = 2
        else
            out[i,j] = 3
        end
    end
    return out
end
function get_color_matrix(grid::GridSpace)
    type_mat = get_type_matrix(grid)
    max_pop  = Float32(max(1, maximum(grid.total_pop)))
    out      = fill(RGBA{Float32}(1f0, 1f0, 1f0, 0f0), grid.rows, grid.cols)

    @inbounds for j in 1:grid.cols, i in 1:grid.rows
        type    = type_mat[i,j]
        opacity = Float32(grid.total_pop[i,j]) / max_pop

        out[i,j] = if type == 0
            RGBA{Float32}(1f0, 1f0, 1f0, 0f0)
        elseif type == 1  # blue
            RGBA{Float32}(0f0, 0f0, 1f0, opacity)
        elseif type == 2  # green
            RGBA{Float32}(0f0, 1f0, 0f0, opacity)
        elseif type == 3  # orange at high density, yellow at low density
            RGBA{Float32}(1f0, 1f0 - 0.5f0*opacity, 0f0, opacity)
        else              # grey (sensitive)
            RGBA{Float32}(0.3f0, 0.3f0, 0.3f0, opacity)
        end
    end
    return out
end

function initialize_grid!(grid::GridSpace;
                          population_fractions::Vector{Float64} = [0.2, 0.2, 0.2],
                          population_types::Vector{Vector{Int}} = [
                              [3000, 3000, 0, 0, 0, 0, 0, 0],
                              [3000, 0, 3000, 0, 0, 0, 0, 0],
                              [0, 0, 0, 0, 3000, 3000, 0, 0]
                          ],
                          seed::Int = 0)

    rng = seed == 0 ? Random.default_rng() : MersenneTwister(seed)
    n_patches   = grid.rows * grid.cols
    assignments = zeros(Int, n_patches)
    current_idx = 1

    for (type_idx, frac) in enumerate(population_fractions)
        n_this = round(Int, frac * n_patches)
        for _ in 1:n_this
            current_idx > n_patches && break
            assignments[current_idx] = type_idx
            current_idx += 1
        end
    end

    Random.shuffle!(rng, assignments)

    @inbounds for idx in 1:n_patches
        i, j = divrem(idx - 1, grid.cols) .+ (1, 1)
        a = assignments[idx]
        a == 0 && continue
        pop = population_types[a]
        grid.nS0[i,j]  = pop[1]; grid.nSR1[i,j] = pop[2]
        grid.nSR2[i,j] = pop[3]; grid.nSS[i,j]  = pop[4]
        grid.nR0[i,j]  = pop[5]; grid.nRS[i,j]  = pop[6]
        grid.nRR1[i,j] = pop[7]; grid.nRR2[i,j] = pop[8]
    end

    update_totals!(grid)
end

## model params
struct ModelParameters
    λ::Float64      # base growth rate
    λP::Float64     # growth rate with sensitive plasmid (nSS, nRS)
    λR::Float64     # growth rate with chrom. resistance, no plasmid (nR0)
    λRP::Float64    # growth rate with resistance + plasmid (nSR1, nSR2)
    λRRP::Float64   # growth rate with chrom. resistance + plasmid (nRS, nRR1, nRR2); no double cost

    s_ss::Float64   # segregation rate, sensitive plasmid
    β_ss::Float64   # conjugation rate, sensitive plasmid
    s1::Float64; β1::Float64  # resistant plasmid 1
    s2::Float64; β2::Float64  # resistant plasmid 2

    α::Float64
    dose::Float64
    K::Float64

    migration_rate::Float64
    extinction_rate::Float64

    dt::Float64   # Gillespie window per step

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

function init_ab_phase(rows::Int, cols::Int, params::ModelParameters; seed::Int = 0)
    A = fill(Int8(0), rows, cols)
    if params.out_of_phase && params.phase_fraction > 0
        rng = MersenneTwister(seed == 0 ? rand(UInt) : UInt(seed))
        @inbounds for j in 1:cols, i in 1:rows
            A[i,j] = rand(rng) < params.phase_fraction ? Int8(1) : Int8(0)
        end
    end
    return A
end

function GridDynamics(grid::GridSpace, params::ModelParameters; seed::Int = 1234)
    rows, cols = grid.rows, grid.cols
    rngs    = [MersenneTwister(seed + i) for i in 1:max(1, nthreads())]
    phase   = init_ab_phase(rows, cols, params; seed = seed + 1000)
    mig_rng = MersenneTwister(seed + 2000)
    states  = [Vector{Int32}(undef, 8) for _ in 1:nthreads()]
    return GridDynamics(grid, params, 0.0, phase, rngs, mig_rng, states, TimeSeriesData())
end

@inline function get_A_patch(t::Float64, phase::Int8, dose::Float64,
                              AB_period::Float64, no_AB_period::Float64)
    AB_period == 0.0 && return 0.0
    cycle = AB_period + no_AB_period
    pos   = mod(t, cycle)
    phase == 0 ? (pos >= no_AB_period ? dose : 0.0) :
                 (pos < AB_period     ? dose : 0.0)
end

@inline function time_to_next_ab_switch(t::Float64, phase::Int8,
                                        AB_period::Float64, no_AB_period::Float64)
    AB_period == 0.0 && return Inf
    cycle = AB_period + no_AB_period
    pos   = mod(t, cycle)
    phase == 0 ? (pos < no_AB_period ? no_AB_period - pos : cycle - pos) :
                 (pos < AB_period    ? AB_period - pos    : cycle - pos)
end

## gillespie
function gillespie!(x::Vector{Int32}, params::ModelParameters, tf::Float64,
                    model_time::Float64, ab_phase::Int8, rng::AbstractRNG)

    λ, λP, λR, λRP, λRRP = params.λ, params.λP, params.λR, params.λRP, params.λRRP
    s_ss, β_ss = params.s_ss, params.β_ss
    s1, β1     = params.s1, params.β1
    s2, β2     = params.s2, params.β2
    α, K, dose   = params.α, params.K, params.dose
    AB_period    = params.AB_period
    no_AB_period = params.no_AB_period

    nS0, nSR1, nSR2, nSS, nR0, nRS, nRR1, nRR2 = x

    invK   = 1.0 / K
    β_ss_K = β_ss * invK
    β1_K   = β1   * invK
    β2_K   = β2   * invK

    dt = 0.0

    while dt < tf
        N   = nS0 + nSR1 + nSR2 + nSS + nR0 + nRS + nRR1 + nRR2
        A   = get_A_patch(model_time + dt, ab_phase, dose, AB_period, no_AB_period)
        αNK = α * N * invK

        donors_r1 = nSR1 + nRR1
        donors_r2 = nSR2 + nRR2
        donors_ss = nSS  + nRS

        total_rate = (
            λ*nS0 + (αNK + A)*nS0 +
            β1_K*nS0*donors_r1 + β2_K*nS0*donors_r2 + β_ss_K*nS0*donors_ss +
            λRP*nSR1 + αNK*nSR1 + s1*nSR1*λRP +
            λRP*nSR2 + αNK*nSR2 + s2*nSR2*λRP +
            λP*nSS + (αNK + A)*nSS + s_ss*nSS*λP +
            λR*nR0 + αNK*nR0 +
            β1_K*nR0*donors_r1 + β2_K*nR0*donors_r2 + β_ss_K*nR0*donors_ss +
            λRRP*nRS + αNK*nRS + s_ss*nRS*λRRP +
            λRRP*nRR1 + αNK*nRR1 + s1*nRR1*λRRP +
            λRRP*nRR2 + αNK*nRR2 + s2*nRR2*λRRP
        )

        total_rate <= 0.0 && break

        τ = -log(rand(rng)) / total_rate

        t_switch = time_to_next_ab_switch(model_time + dt, ab_phase, AB_period, no_AB_period)
        if τ > t_switch && (dt + t_switch) < tf
            dt += t_switch
            continue
        end

        dt += τ
        dt > tf && break

        ru  = rand(rng) * total_rate
        cum = 0.0

        cum += λ*nS0;                ru <= cum && (nS0 += 1; continue)
        cum += (αNK + A)*nS0;        ru <= cum && (nS0 = max(0, nS0-1); continue)
        cum += β1_K*nS0*donors_r1;   ru <= cum && (nS0 = max(0, nS0-1); nSR1 += 1; continue)
        cum += β2_K*nS0*donors_r2;   ru <= cum && (nS0 = max(0, nS0-1); nSR2 += 1; continue)
        cum += β_ss_K*nS0*donors_ss; ru <= cum && (nS0 = max(0, nS0-1); nSS  += 1; continue)

        cum += λRP*nSR1;             ru <= cum && (nSR1 += 1; continue)
        cum += αNK*nSR1;             ru <= cum && (nSR1 = max(0, nSR1-1); continue)
        cum += s1*nSR1*λRP;          ru <= cum && (nSR1 = max(0, nSR1-1); nS0 += 1; continue)

        cum += λRP*nSR2;             ru <= cum && (nSR2 += 1; continue)
        cum += αNK*nSR2;             ru <= cum && (nSR2 = max(0, nSR2-1); continue)
        cum += s2*nSR2*λRP;          ru <= cum && (nSR2 = max(0, nSR2-1); nS0 += 1; continue)

        cum += λP*nSS;               ru <= cum && (nSS += 1; continue)
        cum += (αNK + A)*nSS;        ru <= cum && (nSS = max(0, nSS-1); continue)
        cum += s_ss*nSS*λP;          ru <= cum && (nSS = max(0, nSS-1); nS0 += 1; continue)

        cum += λR*nR0;               ru <= cum && (nR0 += 1; continue)
        cum += αNK*nR0;              ru <= cum && (nR0 = max(0, nR0-1); continue)
        cum += β1_K*nR0*donors_r1;   ru <= cum && (nR0 = max(0, nR0-1); nRR1 += 1; continue)
        cum += β2_K*nR0*donors_r2;   ru <= cum && (nR0 = max(0, nR0-1); nRR2 += 1; continue)
        cum += β_ss_K*nR0*donors_ss; ru <= cum && (nR0 = max(0, nR0-1); nRS  += 1; continue)

        cum += λRRP*nRS;             ru <= cum && (nRS += 1; continue)
        cum += αNK*nRS;              ru <= cum && (nRS = max(0, nRS-1); continue)
        cum += s_ss*nRS*λRRP;        ru <= cum && (nRS = max(0, nRS-1); nR0 += 1; continue)

        cum += λRRP*nRR1;            ru <= cum && (nRR1 += 1; continue)
        cum += αNK*nRR1;             ru <= cum && (nRR1 = max(0, nRR1-1); continue)
        cum += s1*nRR1*λRRP;         ru <= cum && (nRR1 = max(0, nRR1-1); nR0 += 1; continue)

        cum += λRRP*nRR2;            ru <= cum && (nRR2 += 1; continue)
        cum += αNK*nRR2;             ru <= cum && (nRR2 = max(0, nRR2-1); continue)
        nRR2 = max(0, nRR2-1); nR0 += 1  # final event
    end

    x[1]=nS0; x[2]=nSR1; x[3]=nSR2; x[4]=nSS
    x[5]=nR0; x[6]=nRS;  x[7]=nRR1; x[8]=nRR2
    return x
end

function step!(sim::GridDynamics)
    g, p      = sim.grid, sim.params
    rows, cols = g.rows, g.cols
    n_patches  = rows * cols

    @threads :static for patch_idx in 1:n_patches
        i = div(patch_idx - 1, cols) + 1
        j = mod(patch_idx - 1, cols) + 1
        g.total_pop[i,j] == 0 && continue

        rng = sim.rngs[threadid()]
        buf = sim.thread_states[threadid()]

        @inbounds begin
            buf[1] = g.nS0[i,j];  buf[2] = g.nSR1[i,j]
            buf[3] = g.nSR2[i,j]; buf[4] = g.nSS[i,j]
            buf[5] = g.nR0[i,j];  buf[6] = g.nRS[i,j]
            buf[7] = g.nRR1[i,j]; buf[8] = g.nRR2[i,j]
        end

        gillespie!(buf, p, p.dt, sim.model_time, sim.ab_phase[i,j], rng)

        @inbounds begin
            g.nS0[i,j]  = max(0, buf[1]); g.nSR1[i,j] = max(0, buf[2])
            g.nSR2[i,j] = max(0, buf[3]); g.nSS[i,j]  = max(0, buf[4])
            g.nR0[i,j]  = max(0, buf[5]); g.nRS[i,j]  = max(0, buf[6])
            g.nRR1[i,j] = max(0, buf[7]); g.nRR2[i,j] = max(0, buf[8])
        end
    end

    update_totals!(g)
    handle_extinctions!(sim)
    handle_migration!(sim)
    sim.model_time += p.dt
    return sim
end

@inline function handle_extinctions!(sim::GridDynamics)
    p_ext = sim.params.extinction_rate * sim.params.dt
    p_ext <= 0.0 && return
    g = sim.grid; rng = sim.mig_rng
    @inbounds for j in 1:g.cols, i in 1:g.rows
        if g.total_pop[i,j] > 0 && rand(rng) < p_ext
            g.nS0[i,j] = 0; g.nSR1[i,j] = 0; g.nSR2[i,j] = 0; g.nSS[i,j]  = 0
            g.nR0[i,j] = 0; g.nRS[i,j]  = 0; g.nRR1[i,j] = 0; g.nRR2[i,j] = 0
            g.total_pop[i,j] = 0
        end
    end
end

@inline lin_idx(i, j, cols) = (i - 1) * cols + j
@inline cart_idx(idx, cols) = (div(idx - 1, cols) + 1, mod(idx - 1, cols) + 1)

function handle_migration!(sim::GridDynamics)
    r, c = sim.grid.rows, sim.grid.cols
    g, p = sim.grid, sim.params
    p.migration_rate <= 0.0 && return

    rng       = sim.mig_rng
    events    = Vector{NTuple{10,Int32}}()
    dir_probs = [0.25, 0.25, 0.25, 0.25]
    sizehint!(events, r * c * 2)

    @inbounds for j in 1:c, i in 1:r
        g.total_pop[i,j] == 0 && continue

        nbrs = ((i, mod1(j+1, c)), (i, mod1(j-1, c)),
                (mod1(i+1, r), j), (mod1(i-1, r), j))

        mig  = p.migration_rate * p.dt

        nS0  = g.nS0[i,j];  nSR1 = g.nSR1[i,j]; nSR2 = g.nSR2[i,j]; nSS  = g.nSS[i,j]
        nR0  = g.nR0[i,j];  nRS  = g.nRS[i,j];  nRR1 = g.nRR1[i,j]; nRR2 = g.nRR2[i,j]

        mS0_tot  = nS0  > 0 ? min(rand(rng, Poisson(nS0  * mig)), nS0)  : 0
        mSR1_tot = nSR1 > 0 ? min(rand(rng, Poisson(nSR1 * mig)), nSR1) : 0
        mSR2_tot = nSR2 > 0 ? min(rand(rng, Poisson(nSR2 * mig)), nSR2) : 0
        mSS_tot  = nSS  > 0 ? min(rand(rng, Poisson(nSS  * mig)), nSS)  : 0
        mR0_tot  = nR0  > 0 ? min(rand(rng, Poisson(nR0  * mig)), nR0)  : 0
        mRS_tot  = nRS  > 0 ? min(rand(rng, Poisson(nRS  * mig)), nRS)  : 0
        mRR1_tot = nRR1 > 0 ? min(rand(rng, Poisson(nRR1 * mig)), nRR1) : 0
        mRR2_tot = nRR2 > 0 ? min(rand(rng, Poisson(nRR2 * mig)), nRR2) : 0

        mS0_d  = mS0_tot  > 0 ? rand(rng, Multinomial(mS0_tot,  dir_probs)) : [0,0,0,0]
        mSR1_d = mSR1_tot > 0 ? rand(rng, Multinomial(mSR1_tot, dir_probs)) : [0,0,0,0]
        mSR2_d = mSR2_tot > 0 ? rand(rng, Multinomial(mSR2_tot, dir_probs)) : [0,0,0,0]
        mSS_d  = mSS_tot  > 0 ? rand(rng, Multinomial(mSS_tot,  dir_probs)) : [0,0,0,0]
        mR0_d  = mR0_tot  > 0 ? rand(rng, Multinomial(mR0_tot,  dir_probs)) : [0,0,0,0]
        mRS_d  = mRS_tot  > 0 ? rand(rng, Multinomial(mRS_tot,  dir_probs)) : [0,0,0,0]
        mRR1_d = mRR1_tot > 0 ? rand(rng, Multinomial(mRR1_tot, dir_probs)) : [0,0,0,0]
        mRR2_d = mRR2_tot > 0 ? rand(rng, Multinomial(mRR2_tot, dir_probs)) : [0,0,0,0]

        for (d, (ni, nj)) in enumerate(nbrs)
            tot = mS0_d[d] + mSR1_d[d] + mSR2_d[d] + mSS_d[d] +
                  mR0_d[d] + mRS_d[d]  + mRR1_d[d] + mRR2_d[d]
            tot > 0 && push!(events,
                (Int32(lin_idx(i,j,c)), Int32(lin_idx(ni,nj,c)),
                 Int32(mS0_d[d]), Int32(mSR1_d[d]), Int32(mSR2_d[d]), Int32(mSS_d[d]),
                 Int32(mR0_d[d]), Int32(mRS_d[d]),  Int32(mRR1_d[d]), Int32(mRR2_d[d])))
        end
    end

    @inbounds for ev in events
        si, sj = cart_idx(Int(ev[1]), c); ti, tj = cart_idx(Int(ev[2]), c)
        m_S0  = min(Int(ev[3]),  g.nS0[si,sj]);  g.nS0[si,sj]  -= m_S0;  g.nS0[ti,tj]  += m_S0
        m_SR1 = min(Int(ev[4]),  g.nSR1[si,sj]); g.nSR1[si,sj] -= m_SR1; g.nSR1[ti,tj] += m_SR1
        m_SR2 = min(Int(ev[5]),  g.nSR2[si,sj]); g.nSR2[si,sj] -= m_SR2; g.nSR2[ti,tj] += m_SR2
        m_SS  = min(Int(ev[6]),  g.nSS[si,sj]);  g.nSS[si,sj]  -= m_SS;  g.nSS[ti,tj]  += m_SS
        m_R0  = min(Int(ev[7]),  g.nR0[si,sj]);  g.nR0[si,sj]  -= m_R0;  g.nR0[ti,tj]  += m_R0
        m_RS  = min(Int(ev[8]),  g.nRS[si,sj]);  g.nRS[si,sj]  -= m_RS;  g.nRS[ti,tj]  += m_RS
        m_RR1 = min(Int(ev[9]),  g.nRR1[si,sj]); g.nRR1[si,sj] -= m_RR1; g.nRR1[ti,tj] += m_RR1
        m_RR2 = min(Int(ev[10]), g.nRR2[si,sj]); g.nRR2[si,sj] -= m_RR2; g.nRR2[ti,tj] += m_RR2
    end

    update_totals!(g)
end

## main sim loop
function run_simulation!(sim::GridDynamics; n_steps::Int = 100)

    update_time_series!(sim.time_series, sim.grid, sim.model_time)

    for step_n in 1:n_steps
        step!(sim)
        update_time_series!(sim.time_series, sim.grid, sim.model_time)

        if is_extinct(sim.grid)
            println("Extinction at step $step_n (t = $(sim.model_time))")
            break
        end
    end

    return sim
end