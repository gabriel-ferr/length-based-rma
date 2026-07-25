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
# --------------------------------------------------------------------------------
results = Array{Float64, 4}(undef, NUM_INITIAL_CONDITIONS, length(TIMESERIES_LENGTH_RANGE), 3, 2)
# --------------------------------------------------------------------------------
@showprogress for i ∈ 1:NUM_INITIAL_CONDITIONS
    for j ∈ eachindex(TIMESERIES_LENGTH_RANGE)
        sys = CoupledODEs(harmonic_oscilator!, rand(2), ho_params)
        ho, _ = trajectory(sys, TIMESERIES_LENGTH_RANGE[j] * 0.01 + 1; Ttr = TRANSIENT_TIME * 0.01, Δt = 0.01)
        ho = ho[1:TIMESERIES_LENGTH_RANGE[j]]

        dmax = sqrt((maximum(pairwise(Euclidean(), ho)) - minimum(pairwise(Euclidean(), ho))).^2)
        # --------------------------------------------------------------------------------
        rp = RecurrenceMatrix(ho, dmax * FIXED_HO_THRESHOLD)
        results[i, j, 1, 1] = dl_average(rp; lmin = 2)
        results[i, j, 2, 1] = trappingtime(rp; lmin = 2)
        results[i, j, 3, 1] = meanrecurrencetime(rp; lmin = 2)

        results[i, j, 1, 2] = complexity(RecurrenceAverageDiagonal(2, dmax * FIXED_HO_THRESHOLD), ho)
        results[i, j, 2, 2] = complexity(RecurrenceTrappingTime(2, dmax * FIXED_HO_THRESHOLD), ho)
        results[i, j, 3, 2] = complexity(RecurrenceMeanTime(2, dmax * FIXED_HO_THRESHOLD), ho)
        # --------------------------------------------------------------------------------
    end
end
# --------------------------------------------------------------------------------
BSON.@save datadir("results_ho.bson") results