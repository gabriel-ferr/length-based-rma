# --------------------------------------------------------------------------------
using DrWatson
@quickactivate "Length based RMA"
include(projectdir("settings.jl"))
# --------------------------------------------------------------------------------
using BenchmarkTools
using BSON
using Distributions
using ProgressMeter
using RecurrenceAnalysis
using RecurrenceMicrostatesAnalysis
# --------------------------------------------------------------------------------
include(srcdir("measures.jl"))
# --------------------------------------------------------------------------------
benchmarks_rp_with_d = []
benchmarks_rp_without_d = []
benchmarks_method_with_d = []
benchmarks_method_without_d = []

benchmarks_rp_with_tt = []
benchmarks_rp_without_tt = []
benchmarks_method_with_tt = []
benchmarks_method_without_tt = []

benchmarks_rp_with_mrt = []
benchmarks_rp_without_mrt = []
benchmarks_method_with_mrt = []
benchmarks_method_without_mrt = []
# --------------------------------------------------------------------------------
function compute_from_rp_with_d(data)
    rp = RecurrenceMatrix(data, FIXED_THRESHOLD)
    return dl_average(rp.data; lmin = 2, theiler=0)
end

function compute_from_rp_without_d(rp)
    return dl_average(rp.data; lmin = 2, theiler=0)
end

function compute_from_rp_with_tt(data)
    rp = RecurrenceMatrix(data, FIXED_THRESHOLD)
    return trappingtime(rp; lmin = 2, theiler=0)
end

function compute_from_rp_without_tt(rp)
    return trappingtime(rp; lmin = 2, theiler=0)
end

function compute_from_rp_with_mrt(data)
    rp = RecurrenceMatrix(data, FIXED_THRESHOLD)
    return mean_recurrence_time_v2(rp; wmin = 2)
end

function compute_from_rp_without_mrt(rp)
    return mean_recurrence_time_v2(rp; wmin = 2)
end

# --------------------------------------------------------------------------------

function compute_from_probs_with_d(data)
    return complexity(RecurrenceAverageDiagonal(2, FIXED_THRESHOLD), data)
end

function compute_from_probs_without_d(probs, rmspace)
    return RecurrenceMicrostatesAnalysis.measure(RecurrenceAverageDiagonal(2, FIXED_THRESHOLD), rmspace, probs)
end

function compute_from_probs_with_tt(data)
    return complexity(RecurrenceTrappingTime(2, FIXED_THRESHOLD), data)
end

function compute_from_probs_without_tt(probs, rmspace)
    return RecurrenceMicrostatesAnalysis.measure(RecurrenceTrappingTime(2, FIXED_THRESHOLD), rmspace, probs)
end

function compute_from_probs_with_mrt(data)
    return complexity(RecurrenceMeanTime(2, FIXED_THRESHOLD), data)
end

function compute_from_probs_without_mrt(probs, rmspace)
    return RecurrenceMicrostatesAnalysis.measure(RecurrenceMeanTime(2, FIXED_THRESHOLD), rmspace, probs)
end
# --------------------------------------------------------------------------------
base_rmspace = RecurrenceMicrostates(FIXED_THRESHOLD, 3)
# --------------------------------------------------------------------------------
@showprogress for i ∈ eachindex(TIMESERIES_LENGTH_RANGE)
    global input = StateSpaceSet(rand(Uniform(0, 1), TIMESERIES_LENGTH_RANGE[i]))
    global base_rp = RecurrenceMatrix(input, FIXED_THRESHOLD)
    global base_probs = probabilities(base_rmspace, input)

    push!(benchmarks_rp_with_d, @benchmark compute_from_rp_with_d(input) seconds = 10 samples = 100)
    push!(benchmarks_rp_without_d, @benchmark compute_from_rp_without_d(base_rp) seconds = 10 samples = 100)
    push!(benchmarks_method_with_d, @benchmark compute_from_probs_with_d(input) seconds = 10 samples = 100)
    push!(benchmarks_method_without_d, @benchmark compute_from_probs_without_d(base_probs, base_rmspace) seconds = 10 samples = 100)
    push!(benchmarks_rp_with_tt, @benchmark compute_from_rp_with_tt(input) seconds = 10 samples = 100)
    push!(benchmarks_rp_without_tt, @benchmark compute_from_rp_without_tt(base_rp) seconds = 10 samples = 100)
    push!(benchmarks_method_with_tt, @benchmark compute_from_probs_with_tt(input) seconds = 10 samples = 100)
    push!(benchmarks_method_without_tt, @benchmark compute_from_probs_without_tt(base_probs, base_rmspace) seconds = 10 samples = 100)
    push!(benchmarks_rp_with_mrt, @benchmark compute_from_rp_with_mrt(input) seconds = 10 samples = 100)
    push!(benchmarks_rp_without_mrt, @benchmark compute_from_rp_without_mrt(base_rp) seconds = 10 samples = 100)
    push!(benchmarks_method_with_mrt, @benchmark compute_from_probs_with_mrt(input) seconds = 10 samples = 100)
    push!(benchmarks_method_without_mrt, @benchmark compute_from_probs_without_mrt(base_probs, base_rmspace) seconds = 10 samples = 100)
end
# --------------------------------------------------------------------------------
BSON.@save datadir("performance.bson") benchmarks_rp_with_d benchmarks_rp_without_d benchmarks_method_with_d benchmarks_method_without_d benchmarks_rp_with_tt benchmarks_rp_without_tt benchmarks_method_with_tt benchmarks_method_without_tt benchmarks_rp_with_mrt benchmarks_rp_without_mrt benchmarks_method_with_mrt benchmarks_method_without_mrt
