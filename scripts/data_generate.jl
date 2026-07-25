# --------------------------------------------------------------------------------
#       This script generates the time-series for the systems used.
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
using Distributions
using DynamicalSystemsBase
using PredefinedDynamicalSystems
using ProgressMeter
# --------------------------------------------------------------------------------
#   Gaussian white noise (GWN)
data = Matrix{StateSpaceSet}(undef, NUM_INITIAL_CONDITIONS, length(TIMESERIES_LENGTH_RANGE))
@showprogress for j ∈ eachindex(TIMESERIES_LENGTH_RANGE)
    for i ∈ 1:NUM_INITIAL_CONDITIONS
        data[i, j] = StateSpaceSet(rand(Normal(0, 1), TIMESERIES_LENGTH_RANGE[j]))
    end
end

BSON.@save datadir("raw", "gwn.bson") data
# --------------------------------------------------------------------------------
#   Red noise
data = Matrix{StateSpaceSet}(undef, NUM_INITIAL_CONDITIONS, length(TIMESERIES_LENGTH_RANGE))
rng = Oof2RNG(GaussRNG(), 1e-3, 1.0, 1e2)
@showprogress for j ∈ eachindex(TIMESERIES_LENGTH_RANGE)
    for i ∈ 1:NUM_INITIAL_CONDITIONS
        data[i, j] = StateSpaceSet([randoof2(rng) for i in 1:TIMESERIES_LENGTH_RANGE[j]])
    end
end

BSON.@save datadir("raw", "red_noise.bson") data
# --------------------------------------------------------------------------------
#   Hénon map
data = Matrix{StateSpaceSet}(undef, NUM_INITIAL_CONDITIONS, length(TIMESERIES_LENGTH_RANGE))
@showprogress for j ∈ eachindex(TIMESERIES_LENGTH_RANGE)
    for i ∈ 1:NUM_INITIAL_CONDITIONS
        while true
            try 
                sys = PredefinedDynamicalSystems.henon(rand(Uniform(0, 1), 2); a = henon_params[1], b = henon_params[2])
                R, _ = trajectory(sys, TIMESERIES_LENGTH_RANGE[j]; Ttr = TRANSIENT_TIME)
                data[i, j] = R[1:TIMESERIES_LENGTH_RANGE[j]]
                break
            catch
                continue
            end
        end
    end
end
BSON.@save datadir("raw", "henon.bson") data
# --------------------------------------------------------------------------------
#   Lorenz system
data = Matrix{StateSpaceSet}(undef, NUM_INITIAL_CONDITIONS, length(TIMESERIES_LENGTH_RANGE))
@showprogress for j ∈ eachindex(TIMESERIES_LENGTH_RANGE)
    for i ∈ 1:NUM_INITIAL_CONDITIONS
        sys = PredefinedDynamicalSystems.lorenz(rand(Uniform(0, 1), 3); σ = lorenz_params[1], β = lorenz_params[2], ρ = lorenz_params[3])
        R, _ = trajectory(sys, TIMESERIES_LENGTH_RANGE[j] * 0.2 + 1; Ttr = TRANSIENT_TIME * 0.2, Δt = 0.2)
        data[i, j] = R[1:TIMESERIES_LENGTH_RANGE[j]]
    end
end
BSON.@save datadir("raw", "lorenz.bson") data
# --------------------------------------------------------------------------------
#   Harmonic oscillator
function harmonic_oscilator!(du, u, p, _)
    ω = p
    x, v = u
    du[1] = v
    du[2] = -(ω^2) * x
end

data = Matrix{StateSpaceSet}(undef, NUM_INITIAL_CONDITIONS, length(TIMESERIES_LENGTH_RANGE))
@showprogress for j ∈ eachindex(TIMESERIES_LENGTH_RANGE)
    for i ∈ 1:NUM_INITIAL_CONDITIONS
        sys = CoupledODEs(harmonic_oscilator!, rand(2), ho_params)
        R, _ = trajectory(sys, TIMESERIES_LENGTH_RANGE[j] * 0.01 + 1; Ttr = TRANSIENT_TIME * 0.01, Δt = 0.01)
        data[i, j] = StateSpaceSet(R[1:TIMESERIES_LENGTH_RANGE[j], 1])
    end
end
BSON.@save datadir("raw", "harmonic_oscilator.bson") data