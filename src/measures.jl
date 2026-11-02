#
#       This code extends the package RecurrenceMicrostatesAnalysis.jl
#   by implementing the proposed estimators.
# =============================================================================
using ComplexityMeasures
using RecurrenceMicrostatesAnalysis
# =============================================================================
#       Structures
struct RecurrenceAverageDiagonal{L, M, SM} <: ComplexityEstimator 
    ε::Real
    metric::M
    sampling::SM
end

struct RecurrenceMeanTime{L, M, SM} <: ComplexityEstimator 
    ε::Real
    metric::M
    sampling::SM
end

struct RecurrenceTrappingTime{L, M, SM} <: ComplexityEstimator 
    ε::Real
    metric::M
    sampling::SM
end

# -- Constructors
RecurrenceAverageDiagonal(Lmin::Int, ε::Real; metric::M = RecurrenceMicrostatesAnalysis.DEFAULT_METRIC, ratio::Float64 = 0.1, sampling::SM = SRandom(ratio)) where {M,SM} = RecurrenceAverageDiagonal{Lmin,M,SM}(ε, metric, sampling)
RecurrenceMeanTime(Lmin::Int, ε::Real; metric::M = RecurrenceMicrostatesAnalysis.DEFAULT_METRIC, ratio::Float64 = 0.1, sampling::SM = SRandom(ratio)) where {M,SM} = RecurrenceMeanTime{Lmin,M,SM}(ε, metric, sampling)
RecurrenceTrappingTime(Lmin::Int, ε::Real; metric::M = RecurrenceMicrostatesAnalysis.DEFAULT_METRIC, ratio::Float64 = 0.1, sampling::SM = SRandom(ratio)) where {M,SM} = RecurrenceTrappingTime{Lmin,M,SM}(ε, metric, sampling)
# =============================================================================
#   Complexity functions
function ComplexityMeasures.complexity(
        c::RecurrenceAverageDiagonal{L}, 
        x::Union{StateSpaceSet, Vector{<:Real}}, 
        y::Union{StateSpaceSet, Vector{<:Real}} = x;
    ) where {L}
    #       Assert minimum diagonal length.
    @assert 1 ≤ L ≤ 2 "The minimum diagonal length must be 1 or 2."
    return measure(c, x, y)
end

function ComplexityMeasures.complexity(
        c::RecurrenceMeanTime{L}, 
        x::Union{StateSpaceSet, Vector{<:Real}}, 
        y::Union{StateSpaceSet, Vector{<:Real}} = x;
    ) where {L}
    #       Assert minimum diagonal length.
    @assert 1 ≤ L ≤ 2 "The minimum diagonal length must be 1 or 2."
    return measure(c, x, y)
end

function ComplexityMeasures.complexity(
        c::RecurrenceTrappingTime{L}, 
        x::Union{StateSpaceSet, Vector{<:Real}}, 
        y::Union{StateSpaceSet, Vector{<:Real}} = x;
    ) where {L}
    #       Assert minimum diagonal length.
    @assert 1 ≤ L ≤ 2 "The minimum diagonal length must be 1 or 2."
    return measure(c, x, y)
end
# =============================================================================
#   Measure functions
#   L_min = 1
function measure(
    opt::RecurrenceAverageDiagonal{1}, 
    x::Union{StateSpaceSet, Vector{<:Real}},
    y::Union{StateSpaceSet, Vector{<:Real}} = x;
    )
    @assert length(x) == length(y) "Input must have same size."

    #   Compute the recurrence microstates distributions
    cnts_overall, outs_overall = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, DiagonalMicrostate(2); metric = opt.metric, sampling = opt.sampling
        ),
        x, y
    )
    probs_overall = Probabilities(cnts_overall, outs_overall)

    cnts_last_col, _ = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, 1; metric = opt.metric, sampling = Full(),
            space = LimitedRectSamplingSpace(length(x), 1, 0, length(y) - 1),
        ),
        x, y
    )

    #   Compute quantities
    RR = RecurrenceMicrostatesAnalysis.measure(RecurrenceRate(opt.ε), probs_overall)
    ∂D1 = (2 * cnts_last_col[2]) + 1

    return ((length(x)^2) / ((length(x) - 1)^2)) * (RR) / (probs_overall[2] + (∂D1 / ((length(x) - 1)^2)))
end

function measure(
    opt::RecurrenceTrappingTime{1}, 
    x::Union{StateSpaceSet, Vector{<:Real}},
    y::Union{StateSpaceSet, Vector{<:Real}} = x;
    )
    @assert length(x) == length(y) "Input must have same size."

    #   Compute the recurrence microstates distributions
    cnts_overall, outs_overall = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, RectMicrostate(1, 2); metric = opt.metric, sampling = opt.sampling
        ),
        x, y
    )
    probs_overall = Probabilities(cnts_overall, outs_overall)

    cnts_last_col, _ = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, 1; metric = opt.metric, sampling = Full(),
            space = LimitedRectSamplingSpace(length(x), 1, 0, length(y) - 1),
        ),
        x, y
    )

    #   Compute quantities
    RR = RecurrenceMicrostatesAnalysis.measure(RecurrenceRate(opt.ε), probs_overall)
    ∂TT1 = cnts_last_col[2]

    return ((length(x)^2) / ((length(x) - 1)^2)) * (RR) / (probs_overall[2] + (∂TT1 / ((length(x) - 1)^2)))
end

function measure(
    opt::RecurrenceMeanTime{1}, 
    x::Union{StateSpaceSet, Vector{<:Real}},
    y::Union{StateSpaceSet, Vector{<:Real}} = x;
    )
    @assert length(x) == length(y) "Input must have same size."

    #   Compute the recurrence microstates distributions
    cnts_overall, outs_overall = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, RectMicrostate(1, 2); metric = opt.metric, sampling = opt.sampling
        ),
        x, y
    )
    probs_overall = Probabilities(cnts_overall, outs_overall)

    cnts_last_col, _ = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, 1; metric = opt.metric, sampling = Full(),
            space = LimitedRectSamplingSpace(length(x), 1, 0, length(y) - 1),
        ),
        x, y
    )

    #   Compute quantities
    RR = RecurrenceMicrostatesAnalysis.measure(RecurrenceRate(opt.ε), probs_overall)
    ∂MT1 = cnts_last_col[1]

    return ((length(x)^2) / ((length(x) - 1)^2)) * (1 - RR) / (probs_overall[3] + (∂MT1 / ((length(x) - 1)^2)))
end

#################
#   L_min = 2
function measure(
    opt::RecurrenceAverageDiagonal{2}, 
    x::Union{StateSpaceSet, Vector{<:Real}},
    y::Union{StateSpaceSet, Vector{<:Real}} = x;
    )
    @assert length(x) == length(y) "Input must have same size."

    #   Compute the recurrence microstates distributions
    cnts_overall, outs_overall = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, DiagonalMicrostate(3); metric = opt.metric, sampling = opt.sampling
        ),
        x, y
    )
    probs_overall = Probabilities(cnts_overall, outs_overall)

    cnts_last_col, _ = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, DiagonalMicrostate(2); metric = opt.metric, sampling = Full(),
            space = LimitedRectSamplingSpace(1, 1, 0, 1),
        ),
        x, y
    )

    #   Compute quantities
    ∂D2 = (2 * cnts_last_col[4]) - 1

    return 1 + ((probs_overall[4] + probs_overall[8]) / (probs_overall[7] + (∂D2 / ((length(x) - 2)^2))))
end

function measure(
    opt::RecurrenceTrappingTime{2}, 
    x::Union{StateSpaceSet, Vector{<:Real}},
    y::Union{StateSpaceSet, Vector{<:Real}} = x;
    )
    @assert length(x) == length(y) "Input must have same size."

    #   Compute the recurrence microstates distributions
    cnts_overall, outs_overall = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, RectMicrostate(3, 1); metric = opt.metric, sampling = opt.sampling
        ),
        x, y
    )
    probs_overall = Probabilities(cnts_overall, outs_overall)

    cnts_last_col, _ = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, RectMicrostate(2, 1); metric = opt.metric, sampling = Full(),
            space = LimitedRectSamplingSpace(1, 1, 0, 1),
        ),
        x, y
    )

    #   Compute quantities
    ∂TT2 = cnts_last_col[4]

    return 1 + ((probs_overall[4] + probs_overall[8]) / (probs_overall[7] + (∂TT2 / ((length(x) - 2)^2))))
end

function measure(
    opt::RecurrenceMeanTime{2}, 
    x::Union{StateSpaceSet, Vector{<:Real}},
    y::Union{StateSpaceSet, Vector{<:Real}} = x;
    )
    @assert length(x) == length(y) "Input must have same size."

    #   Compute the recurrence microstates distributions
    cnts_overall, outs_overall = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, RectMicrostate(3, 1); metric = opt.metric, sampling = opt.sampling
        ),
        x, y
    )
    probs_overall = Probabilities(cnts_overall, outs_overall)

    cnts_last_col, _ = counts_and_outcomes(
        RecurrenceMicrostates(
            opt.ε, RectMicrostate(2, 1); metric = opt.metric, sampling = Full(),
            space = LimitedRectSamplingSpace(1, 1, 0, 1),
        ),
        x, y
    )

    #   Compute quantities
    ∂MT2 = cnts_last_col[1]

    return 1 + ((probs_overall[1] + probs_overall[5]) / (probs_overall[2] + (∂MT2 / ((length(x) - 2)^2))))
end
# =============================================================================
#   Modify the TriangleSamplingSpace and RectSamplingSpace to select a specific
#   area of the RP. This part is VERY complex, because I need to adapt some backend
#   stuff of the package ... Note that I only implemented it for rectangle microstates,
#   because its implementation is a little hard, and now I just need it for rectangle.
struct LimitedRectSamplingSpace <: AbstractSamplingSpace
    init_x::Int
    init_y::Int
    len_x::Int
    len_y::Int
end

LimitedRectSamplingSpace(x = 1, y = 1, len_x = 0, len_y = 0) = LimitedRectSamplingSpace(x, y, len_x, len_y)
# =============================================================================
struct ExtendedSSRect2 <: RecurrenceMicrostatesAnalysis.InternalSamplingSpace
    init_W::Int
    init_H::Int
    len_W::Int
    len_H::Int
end
# =============================================================================
function RecurrenceMicrostatesAnalysis.SamplingSpace(
        ::RecurrenceMicrostatesAnalysis.Rect2Microstate{W, H},
        space::LimitedRectSamplingSpace,
        x::StateSpaceSet,
        y::StateSpaceSet
    ) where {W, H}

    Lx, Ly = length(x), length(y)

    @assert 1 ≤ space.init_x ≤ Lx - W + 1 "Invalid `init_x`. Use 1 ≤ init_x ≤ $(Lx - W + 1)."
    @assert 1 ≤ space.init_y ≤ Ly - H + 1 "Invalid `init_y`. Use 1 ≤ init_y ≤ $(Ly - H + 1)."

    # número de posições iniciais válidas
    len_W = space.len_x <= 0 ? (Lx - W + 1) - (space.init_x - 1) : space.len_x
    len_H = space.len_y <= 0 ? (Ly - H + 1) - (space.init_y - 1) : space.len_y

    @assert len_W ≥ 1 "Invalid x-window length."
    @assert len_H ≥ 1 "Invalid y-window length."
    @assert space.init_x + len_W + W - 2 ≤ Lx "The x-window exceeds the time-series length. Use len_x ≤ $(Lx - W - space.init_x + 2)."
    @assert space.init_y + len_H + H - 2 ≤ Ly "The y-window exceeds the time-series length. Use len_y ≤ $(Ly - H - space.init_y + 2)."

    return ExtendedSSRect2(space.init_x, space.init_y, len_W, len_H)
end

function RecurrenceMicrostatesAnalysis.SamplingSpace(
        ::DiagonalMicrostate{N},
        space::LimitedRectSamplingSpace,
        x::StateSpaceSet,
        y::StateSpaceSet
    ) where {N}

    Lx, Ly = length(x), length(y)

    @assert 1 ≤ space.init_x ≤ Lx - N + 1 "Invalid `init_x`. Use 1 ≤ init_x ≤ $(Lx - N + 1)."
    @assert 1 ≤ space.init_y ≤ Ly - N + 1 "Invalid `init_y`. Use 1 ≤ init_y ≤ $(Ly - N + 1)."

    # número de posições iniciais válidas
    len_W = space.len_x <= 0 ? (Lx - N + 1) - (space.init_x - 1) : space.len_x
    len_H = space.len_y <= 0 ? (Ly - N + 1) - (space.init_y - 1) : space.len_y

    @assert len_W ≥ 1 "Invalid x-window length."
    @assert len_H ≥ 1 "Invalid y-window length."
    @assert space.init_x + len_W + N - 2 ≤ Lx "The x-window exceeds the time-series length. Use len_x ≤ $(Lx - N - space.init_x + 2)."
    @assert space.init_y + len_H + N - 2 ≤ Ly "The y-window exceeds the time-series length. Use len_y ≤ $(Ly - N - space.init_y + 2)."

    return ExtendedSSRect2(space.init_x, space.init_y, len_W, len_H)
end
# =============================================================================
function RecurrenceMicrostatesAnalysis.get_num_samples(mode::SRandom{<:Real}, space::ExtendedSSRect2)
    return ceil(Int, mode.sampling_factor * space.len_W * space.len_H)
end

function RecurrenceMicrostatesAnalysis.get_num_samples(::Full, space::ExtendedSSRect2)
    return space.len_H * space.len_W
end
# =============================================================================
function RecurrenceMicrostatesAnalysis.get_sample(
        ::RecurrenceMicrostatesAnalysis.CPUCore,
        ::RecurrenceMicrostatesAnalysis.SRandom,
        space::ExtendedSSRect2,
        rng,
        _
    )
    
    i = rand(rng, space.init_W:(space.init_W + space.len_W - 1))
    j = rand(rng, space.init_H:(space.init_H + space.len_H - 1))

    return i, j
end

function RecurrenceMicrostatesAnalysis.get_sample(
        ::RecurrenceMicrostatesAnalysis.CPUCore,
        ::RecurrenceMicrostatesAnalysis.Full,
        space::ExtendedSSRect2,
        _,
        m
)
    i = ((m - 1) % space.len_W) + space.init_W
    j = ((m - 1) ÷ space.len_W) + space.init_H

    return i, j
end
# =============================================================================
# using Distributions
# using Random
# using RecurrenceAnalysis

# include("rmt.jl")

# function test()
#     serie = StateSpaceSet(rand(Uniform(0, 1), 1000))
#     rp = RecurrenceMatrix(serie, 0.98)

#     ref = dl_average(rp; lmin = 1)
#     comp = measure(RecurrenceAverageDiagonal(1, 0.98), serie, serie) 

#     println("D1: ref = $(ref), comp = $(comp); % = $(100 * abs(ref - comp) / ref)")

#     ref = trappingtime(rp; lmin = 1)
#     comp = measure(RecurrenceTrappingTime(1, 0.98), serie, serie)

#     println("TT1: ref = $(ref), comp = $(comp); % = $(100 * abs(ref - comp) / ref)")

#     ref = mean_recurrence_time(rp; wmin = 1)
#     comp = measure(RecurrenceMeanTime(1, 0.98), serie, serie)
#     println("MT1: ref = $(ref), comp = $(comp); % = $(100 * abs(ref - comp) / ref)")

#     ref = dl_average(rp; lmin = 2)
#     comp = measure(RecurrenceAverageDiagonal(2, 0.98), serie, serie) 

#     println("D2: ref = $(ref), comp = $(comp); % = $(100 * abs(ref - comp) / ref)")

#     ref = trappingtime(rp; lmin = 2)
#     comp = measure(RecurrenceTrappingTime(2, 0.98), serie, serie)

#     println("TT2: ref = $(ref), comp = $(comp); % = $(100 * abs(ref - comp) / ref)")

#     ref = mean_recurrence_time(rp; wmin = 2)
#     comp = measure(RecurrenceMeanTime(2, 0.98), serie, serie)
#     println("MT2: ref = $(ref), comp = $(comp); % = $(100 * abs(ref - comp) / ref)")

#     # rng = TaskLocalRNG()
#     # core = RecurrenceMicrostatesAnalysis.CPUCore()
#     # samp = RecurrenceMicrostatesAnalysis.SRandom(0.1)
#     # spc = RecurrenceMicrostatesAnalysis.SamplingSpace(
#     #     DiagonalMicrostate(2),
#     #     LimitedRectSamplingSpace(1, 1, length(serie) - 1, 1),
#     #     serie,
#     #     serie
#     # )

#     # for _ ∈ 1:1000
#     #     i, j = RecurrenceMicrostatesAnalysis.get_sample(core, samp, spc, rng, nothing)
#     #     print("($(i), $(j)); \t")
#     # end
# end

# test()