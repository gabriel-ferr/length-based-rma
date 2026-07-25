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
using CairoMakie
using Colors
using Statistics
using LaTeXStrings
# --------------------------------------------------------------------------------
set_theme!(
    font = "Arial",
    Axis = (
        xgridvisible = false,
        ygridvisible = false,
        topspinevisible = false,
        rightspinevisible = false,
    )
)
# --------------------------------------------------------------------------------
BSON.@load datadir("results_sampling.bson") results
# --------------------------------------------------------------------------------
errors = Array{Float64, 6}(undef, NUM_INITIAL_CONDITIONS, length(THRESHOLD_RANGE_SAMPLING), length(SAMPLING_RATIO_RANGE), 5, 3, 2)
for i ∈ eachindex(THRESHOLD_RANGE_SAMPLING)
    for j ∈ 1:5
        for k ∈ 1:3
            for n ∈ 1:2
                for p ∈ eachindex(SAMPLING_RATIO_RANGE)
                    for m ∈ 1:NUM_INITIAL_CONDITIONS
                        errors[m, i, p, j, k, n] = abs(results[m, i, p, j, k, 2, n] - results[m, i, 1, j, k, 1, n]) / abs(results[m, i, 1, j, k, 1, n])
                        if (isinf(errors[m, i, p, j, k, n])); errors[m, i, p, j, k, n] = NaN; end
                        if (errors[m, i, p, j, k, n] > 1); errors[m, i, p, j, k, n] = NaN; end
                    end
                end
            end
        end
    end
end
# --------------------------------------------------------------------------------
ax_pos = [
    ([1, 1], "\\text{(A) Gaussian white noise}"),
    ([1, 2], "\\text{(B) Red noise}"),
    ([2, 1], "\\text{(C) Hénon map}"),
    ([2, 2], "\\text{(D) Lorenz system}"),
    ([3, 1], "\\text{(E) Harmonic oscillator, } x"),
]
for n ∈ 1:3
    fig = Figure(size = (1240, 1380))
    for s ∈ 1:5
        # ----------------------------------------------------------------------------
        positions = Vector{Int}(undef, 0)
        values = Vector{Float64}(undef, 0)
        dodges = Vector{Int}(undef, 0)

        for i ∈ 1:NUM_INITIAL_CONDITIONS
            for j ∈ eachindex(THRESHOLD_RANGE_SAMPLING)
                for k ∈ eachindex(SAMPLING_RATIO_RANGE)
                    for m ∈ 1:2
                        if (isnan(errors[i, j, k, s, n, m])); continue; end
                        v = log10(errors[i, j, k, s, n, m])

                        if (isfinite(v))
                            push!(positions, k)
                            push!(values, v)
                            push!(dodges, m)
                        end
                    end
                end
            end
        end

        colors = ifelse.(
            dodges .== 1,
            coloralpha(colorant"steelblue3", 0.6),
            coloralpha(colorant"indianred2", 0.6)
        )
        # ----------------------------------------------------------------------------
        ax = Axis(fig[ax_pos[s][1][1], ax_pos[s][1][2]], ylabel = L"\text{Relative error}", xlabel = L"\text{Sampling ratio}", title = LaTeXString(ax_pos[s][2]))
        hlines!(ax, log10(0.05); color = :black, linestyle = :dash)
        violin!(ax, positions, values; dodge = dodges, color = colors, dodge_gap = 0.1)
        boxplot!(ax, positions, values; dodge = dodges, color = (:black, 0.0), strokecolor = :black, strokewidth = 0.8, dodge_gap = 0.1, whiskerwidth = 0.6)

        ylims!(ax, nothing, 0.0)
        ax.xticks = (1:5, string.(SAMPLING_RATIO_RANGE))
        ax.ytickformat = values -> [L"10^{%$(round(v, digits = 1))}" for v in values]

        axislegend(
            ax,
            [
                PolyElement(color = (:steelblue3, 0.6)),
                PolyElement(color = (:indianred2, 0.6)),
                LineElement(color = :black, linestyle = :dash)
            ],
            [
                L"\ell_{\min}=1",
                L"\ell_{\min}=2",
                L"\text{Error}= 5%"
            ];
            position = :lb
        )
        # ----------------------------------------------------------------------------
        # println("n = $n ; s = $s")
    end
    save(plotsdir("sampling_n=$n.png"), fig; px_per_unit = 2)
end