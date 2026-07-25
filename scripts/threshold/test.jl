# --------------------------------------------------------------------------------
#       This script generates NUM_INITIAL_CONDITIONS time-series for each system
#   using DEFAULT_TIMESERIES_LENGTH points per time-series.
#   Here, we generate:
#   - Gaussian white noise (GWN)
#   - Red noise
#   - Hénon map
#   - Lorenz system
#   - Harmonic oscillator
# --------------------------------------------------------------------------------
using DrWatson
@quickactivate "Length based RMA"
include(projectdir("settings.jl"))
# --------------------------------------------------------------------------------
using BSON
using CorrNoise
using Distances
using Distributions
using DynamicalSystemsBase
using PredefinedDynamicalSystems
using ProgressMeter
using RecurrenceAnalysis
using RecurrenceMicrostatesAnalysis
# --------------------------------------------------------------------------------
function harmonic_oscilator!(du, u, p, _)
    ω = p
    x, v = u
    du[1] = v
    du[2] = -(ω^2) * x
end

rng = Oof2RNG(GaussRNG(), 1e-3, 1.0, 1e2)
# --------------------------------------------------------------------------------
results = Array{Float64, 6}(undef, NUM_INITIAL_CONDITIONS, length(THRESHOLD_RANGE), 5, 3, 2, 2)
# --------------------------------------------------------------------------------
@showprogress for i ∈ 1:NUM_INITIAL_CONDITIONS
    gwn = StateSpaceSet(rand(Normal(0, 1), DEFAULT_TIMESERIES_LENGTH))
    rn =  StateSpaceSet([randoof2(rng) for i in 1:DEFAULT_TIMESERIES_LENGTH])

    sys = PredefinedDynamicalSystems.henon(rand(Uniform(0, 1), 2); a = henon_params[1], b = henon_params[2])
    henon, _ = trajectory(sys, DEFAULT_TIMESERIES_LENGTH; Ttr = TRANSIENT_TIME)
    henon = henon[1:DEFAULT_TIMESERIES_LENGTH]

    while true
        if any(isfinite, Iterators.flatten(henon))
            break
        end

        sys = PredefinedDynamicalSystems.henon(rand(Uniform(0, 1), 2); a = henon_params[1], b = henon_params[2])
        henon, _ = trajectory(sys, DEFAULT_TIMESERIES_LENGTH; Ttr = TRANSIENT_TIME)
        henon = henon[1:DEFAULT_TIMESERIES_LENGTH]
    end
    
    sys = PredefinedDynamicalSystems.lorenz(rand(Uniform(0, 1), 3); σ = lorenz_params[1], β = lorenz_params[2], ρ = lorenz_params[3])
    lorenz, _ = trajectory(sys, DEFAULT_TIMESERIES_LENGTH * 0.2 + 1; Ttr = TRANSIENT_TIME * 0.2, Δt = 0.2)
    lorenz = lorenz[1:DEFAULT_TIMESERIES_LENGTH]

    sys = CoupledODEs(harmonic_oscilator!, rand(2), ho_params)
    ho, _ = trajectory(sys, DEFAULT_TIMESERIES_LENGTH * 0.01 + 1; Ttr = TRANSIENT_TIME * 0.01, Δt = 0.01)
    ho = ho[1:DEFAULT_TIMESERIES_LENGTH]

    input = [gwn, rn, henon, lorenz, ho]

    # -----
    dmax = [
        sqrt((maximum(pairwise(Euclidean(), gwn)) - minimum(pairwise(Euclidean(), gwn))).^2),
        sqrt((maximum(pairwise(Euclidean(), rn)) - minimum(pairwise(Euclidean(), rn))).^2),
        sqrt((maximum(pairwise(Euclidean(), henon)) - minimum(pairwise(Euclidean(), henon))).^2),
        sqrt((maximum(pairwise(Euclidean(), lorenz)) - minimum(pairwise(Euclidean(), lorenz))).^2),
        sqrt((maximum(pairwise(Euclidean(), ho)) - minimum(pairwise(Euclidean(), ho))).^2),
    ]
    # -----

    for j ∈ eachindex(THRESHOLD_RANGE)
        for s ∈ eachindex(input)
            rp = RecurrenceMatrix(input[s], dmax[s] * THRESHOLD_RANGE[j])
            # -----
            results[i, j, s, 1, 1, 1] = dl_average(rp; lmin = 1)
            results[i, j, s, 1, 1, 2] = dl_average(rp; lmin = 2)

            results[i, j, s, 2, 1, 1] = trappingtime(rp; lmin = 1)
            results[i, j, s, 2, 1, 2] = trappingtime(rp; lmin = 2)

            results[i, j, s, 3, 1, 1] = meanrecurrencetime(rp; lmin = 1)
            results[i, j, s, 3, 1, 2] = meanrecurrencetime(rp; lmin = 2)
            # -----
            results[i, j, s, 1, 2, 1] = complexity(RecurrenceAverageDiagonal(1, dmax[s] * THRESHOLD_RANGE[j]), input[s])
            results[i, j, s, 1, 2, 2] = complexity(RecurrenceAverageDiagonal(2, dmax[s] * THRESHOLD_RANGE[j]), input[s])

            results[i, j, s, 2, 2, 1] = complexity(RecurrenceTrappingTime(1, dmax[s] * THRESHOLD_RANGE[j]), input[s])
            results[i, j, s, 2, 2, 2] = complexity(RecurrenceTrappingTime(2, dmax[s] * THRESHOLD_RANGE[j]), input[s])

            results[i, j, s, 3, 2, 1] = complexity(RecurrenceMeanTime(1, dmax[s] * THRESHOLD_RANGE[j]), input[s])
            results[i, j, s, 3, 2, 2] = complexity(RecurrenceMeanTime(2, dmax[s] * THRESHOLD_RANGE[j]), input[s])
        end
    end
end
# --------------------------------------------------------------------------------
BSON.@save datadir("results_threshold.bson") results