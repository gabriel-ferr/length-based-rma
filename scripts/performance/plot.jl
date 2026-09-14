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
using ProgressMeter
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
function confidence_interval(x; α = 0.05)
    n = length(x)
    @assert n ≥ 2

    μ = mean(x)
    s = std(x; corrected=true)      # desvio padrão amostral
    t = quantile(TDist(n - 1), 1 - α/2)

    margin = t * s / sqrt(n)

    return (
        μ,
        μ - margin,
        μ + margin,
        margin,
    )
end
# --------------------------------------------------------------------------------
BSON.@load datadir("performance.bson") benchmarks_rp_with_d benchmarks_rp_without_d benchmarks_method_with_d benchmarks_method_without_d benchmarks_rp_with_tt benchmarks_rp_without_tt benchmarks_method_with_tt benchmarks_method_without_tt benchmarks_rp_with_mrt benchmarks_rp_without_mrt benchmarks_method_with_mrt benchmarks_method_without_mrt
input = [
    benchmarks_rp_with_d,           # 1 - D, RP
    benchmarks_rp_with_tt,          # 2 - TT, RP
    benchmarks_rp_with_mrt,         # 3 - MRT, RP
    benchmarks_rp_without_d,        # 4 - D, NO RP
    benchmarks_rp_without_tt,       # 5 - TT, NO RP
    benchmarks_rp_without_mrt,      # 6 - MRT, NO RP
    benchmarks_method_with_d,       # 7 - D, RMA
    benchmarks_method_with_tt,      # 8 - TT, RMA
    benchmarks_method_with_mrt,     # 9 - MRT, RMA
    benchmarks_method_without_d,    # 10 - D, NO PROBS
    benchmarks_method_without_tt,   # 11 - TT, NO PROBS
    benchmarks_method_without_mrt,  # 12 - MRT, NO PROBS
]
# --------------------------------------------------------------------------------
perform = Array{Float64, 3}(undef, 5, length(TIMESERIES_LENGTH_RANGE), length(input))
for i ∈ eachindex(TIMESERIES_LENGTH_RANGE)
    for j ∈ eachindex(input)
        r = confidence_interval(input[j][i].times)
        perform[1, i, j] = r[1] / 1e9
        perform[2, i, j] = r[2] / 1e9
        perform[3, i, j] = r[3] / 1e9
        perform[4, i, j] = input[j][i].allocs / 1e3
        perform[5, i, j] = input[j][i].memory
    end
end
# --------------------------------------------------------------------------------
fig = Figure(size = (1240, 1395))
# --------------------------------------------------------------------------------
ax = Axis(fig[1, 1], xlabel = L"\text{Time series length}", ylabel = L"\text{Computation time}~(s)", title = L"\textbf{(A)}", yscale = log10)

scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[1, :, 1]; marker = :utriangle, markersize = 12, strokecolor = :blue, strokewidth = 1, color = :blue, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[1, :, 2]; marker = :rect, markersize = 12, strokecolor = :orangered1, strokewidth = 1, color = :orangered1, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[1, :, 7]; marker = :utriangle, markersize = 12, strokecolor = :seagreen, strokewidth = 1, color = :seagreen, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[1, :, 8]; marker = :rect, markersize = 12, strokecolor = :violetred4, strokewidth = 1, color = :violetred4, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[1, :, 9]; marker = :circle, markersize = 12, strokecolor = :black, strokewidth = 1, color = :black, markercolor = (:white, 0.0), linewidth = 1, label = L"RT_2,~\text{RMA}")

axislegend(ax, position = :lt)
# --------------------------------------------------------------------------------
ax = Axis(fig[1, 2], xlabel = L"\text{Time series length}", ylabel = L"\text{Computation time}~(s)", title = L"\textbf{(B)}", yscale = log10)

scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[1, :, 4]; marker = :utriangle, markersize = 12, strokecolor = :blue, strokewidth = 1, color = :blue, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[1, :, 5]; marker = :rect, markersize = 12, strokecolor = :orangered1, strokewidth = 1, color = :orangered1, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[1, :, 10]; marker = :utriangle, markersize = 12, strokecolor = :seagreen, strokewidth = 1, color = :seagreen, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[1, :, 11]; marker = :rect, markersize = 12, strokecolor = :violetred4, strokewidth = 1, color = :violetred4, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[1, :, 12]; marker = :circle, markersize = 12, strokecolor = :black, strokewidth = 1, color = :black, markercolor = (:white, 0.0), linewidth = 1, label = L"RT_2,~\text{RMA}")

axislegend(ax, position = :rc)
# --------------------------------------------------------------------------------
ax = Axis(fig[2, 1], xlabel = L"\text{Time series length}", ylabel = L"\text{Memory estimated}", title = L"\textbf{(C)}", yscale = log10, yticks = (
        [2^10, 2^15, 2^20, 2^25, 2^30],
        [L"1\text{KiB}", L"32\text{KiB}", L"1\text{MiB}", L"32\text{MiB}", L"1\text{GiB}"]
    ))
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[5, :, 1]; marker = :utriangle, markersize = 12, strokecolor = :blue, strokewidth = 1, color = :blue, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[5, :, 2]; marker = :rect, markersize = 12, strokecolor = :orangered1, strokewidth = 1, color = :orangered1, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[5, :, 7]; marker = :utriangle, markersize = 12, strokecolor = :seagreen, strokewidth = 1, color = :seagreen, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[5, :, 8]; marker = :rect, markersize = 12, strokecolor = :violetred4, strokewidth = 1, color = :violetred4, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[5, :, 9]; marker = :circle, markersize = 12, strokecolor = :black, strokewidth = 1, color = :black, markercolor = (:white, 0.0), linewidth = 1, label = L"RT_2,~\text{RMA}")

axislegend(ax, position = :rc)
# --------------------------------------------------------------------------------
ax = Axis(fig[2, 2], xlabel = L"\text{Time series length}", ylabel = L"\text{Memory estimated}", title = L"\textbf{(D)}", yscale = log10, yticks = (
        [2^10, 2^15, 2^20, 2^25, 2^30],
        [L"1\text{KiB}", L"32\text{KiB}", L"1\text{MiB}", L"32\text{MiB}", L"1\text{GiB}"]
    ))
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[5, :, 4]; marker = :utriangle, markersize = 12, strokecolor = :blue, strokewidth = 1, color = :blue, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[5, :, 5]; marker = :rect, markersize = 12, strokecolor = :orangered1, strokewidth = 1, color = :orangered1, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[5, :, 10]; marker = :utriangle, markersize = 12, strokecolor = :seagreen, strokewidth = 1, color = :seagreen, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[5, :, 11]; marker = :rect, markersize = 12, strokecolor = :violetred4, strokewidth = 1, color = :violetred4, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[5, :, 12]; marker = :circle, markersize = 12, strokecolor = :black, strokewidth = 1, color = :black, markercolor = (:white, 0.0), linewidth = 1, label = L"RT_2,~\text{RMA}")

axislegend(ax, position = :rc)
# --------------------------------------------------------------------------------
ax = Axis(fig[3, 1], xlabel = L"\text{Time series length}", ylabel = L"\text{Allocations}~(\times 10^3)", title = L"\textbf{(E)}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[4, :, 1]; marker = :utriangle, markersize = 12, strokecolor = :blue, strokewidth = 1, color = :blue, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[4, :, 2]; marker = :rect, markersize = 12, strokecolor = :orangered1, strokewidth = 1, color = :orangered1, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[4, :, 7]; marker = :utriangle, markersize = 12, strokecolor = :seagreen, strokewidth = 1, color = :seagreen, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[4, :, 8]; marker = :rect, markersize = 12, strokecolor = :violetred4, strokewidth = 1, color = :violetred4, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[4, :, 9]; marker = :circle, markersize = 12, strokecolor = :black, strokewidth = 1, color = :black, markercolor = (:white, 0.0), linewidth = 1, label = L"RT_2,~\text{RMA}")

axislegend(ax, position = :lt)
# --------------------------------------------------------------------------------
ax = Axis(fig[3, 2], xlabel = L"\text{Time series length}", ylabel = L"\text{Allocations}~(\times 10^3)", title = L"\textbf{(F)}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[4, :, 4]; marker = :utriangle, markersize = 12, strokecolor = :blue, strokewidth = 1, color = :blue, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[4, :, 5]; marker = :rect, markersize = 12, strokecolor = :orangered1, strokewidth = 1, color = :orangered1, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RP}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[4, :, 10]; marker = :utriangle, markersize = 12, strokecolor = :seagreen, strokewidth = 1, color = :seagreen, markercolor = (:white, 0.0), linewidth = 1, label = L"D_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[4, :, 11]; marker = :rect, markersize = 12, strokecolor = :violetred4, strokewidth = 1, color = :violetred4, markercolor = (:white, 0.0), linewidth = 1, label = L"TT_2,~\text{RMA}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, perform[4, :, 12]; marker = :circle, markersize = 12, strokecolor = :black, strokewidth = 1, color = :black, markercolor = (:white, 0.0), linewidth = 1, label = L"RT_2,~\text{RMA}")

axislegend(ax, position = :lt)
# --------------------------------------------------------------------------------
save(plotsdir("performance.png"), fig; px_per_unit = 6)