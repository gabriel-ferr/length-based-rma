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
using Statistics
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
BSON.@load datadir("results_series_len.bson") results
# --------------------------------------------------------------------------------
errors = Array{Float64, 4}(undef, NUM_INITIAL_CONDITIONS, length(TIMESERIES_LENGTH_RANGE), 4, 3)
condensed = Array{Float64, 4}(undef, length(TIMESERIES_LENGTH_RANGE), 4, 3, 3)

for i ∈ eachindex(TIMESERIES_LENGTH_RANGE)
    for j ∈ 1:3
        for p ∈ 1:4
            for k ∈ 1:NUM_INITIAL_CONDITIONS
                errors[k, i, p, j] = abs(results[k, i, p, j, 2] - results[k, i, p, j, 1]) / abs(results[k, i, p, j, 1])
                if (isinf(errors[k, i, p, j])); errors[k, i, p, j] = NaN; end
                if (errors[k, i, p, j] > 1); errors[k, i, p, j] = NaN; end
            end

            condensed[i, p, j, 1] = mean(errors[:, i, p, j])
            try
                condensed[i, p, j, 2] = quantile(errors[:, i, p, j], 0.1)
                condensed[i, p, j, 3] = quantile(errors[:, i, p, j], 0.9)
            catch
                condensed[i, p, j, 2] = NaN
                condensed[i, p, j, 3] = NaN
            end
        end
    end
end
# --------------------------------------------------------------------------------
fig = Figure(size = (1240, 465))

ax = Axis(fig[1, 1], xlabel = L"\text{Time series length}", ylabel = L"\text{Relative error}", title = L"\textbf{(A)}")
band!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 1, 1, 2], condensed[:, 1, 1, 3]; color = (:blue, 0.1))
band!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 2, 1, 2], condensed[:, 2, 1, 3]; color = (:orangered4, 0.1))
band!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 3, 1, 2], condensed[:, 3, 1, 3]; color = (:seagreen, 0.1))
band!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 4, 1, 2], condensed[:, 4, 1, 3]; color = (:violetred4, 0.1))

scatterlines!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 1, 1, 1]; marker = :utriangle, markersize = 12, strokecolor = :blue, strokewidth = 1, color = :blue, markercolor = (:white, 0.0), label = L"\text{GWN}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 2, 1, 1]; marker = :rect, markersize = 12, strokecolor = :orangered4, strokewidth = 1, color = :orangered4, markercolor = (:white, 0.0), label = L"\text{Hénon map}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 3, 1, 1]; marker = :circle, markersize = 12, strokecolor = :seagreen, strokewidth = 1, color = :seagreen, markercolor = (:white, 0.0), label = L"\text{Lorenz system}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 4, 1, 1]; marker = :cross, markersize = 12, strokecolor = :violetred4, strokewidth = 1, color = :violetred4, markercolor = (:white, 0.0), label = L"\text{Harmonic oscillator, } x")

hlines!(ax, 0.05; color = :black, linestyle = :dash, label = L"\text{Error}=5%")

axislegend(ax, position = :rt)

ax = Axis(fig[1, 2], xlabel = L"\text{Time series length}", ylabel = L"\text{Relative error}", title = L"\textbf{(B)}")
band!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 1, 2, 2], condensed[:, 1, 2, 3]; color = (:blue, 0.1))
band!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 2, 2, 2], condensed[:, 2, 2, 3]; color = (:orangered4, 0.1))
band!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 3, 2, 2], condensed[:, 3, 2, 3]; color = (:seagreen, 0.1))
band!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 4, 2, 2], condensed[:, 4, 2, 3]; color = (:violetred4, 0.1))

scatterlines!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 1, 2, 1]; marker = :utriangle, markersize = 12, strokecolor = :blue, strokewidth = 1, color = :blue, markercolor = (:white, 0.0), label = L"\text{GWN}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 2, 2, 1]; marker = :rect, markersize = 12, strokecolor = :orangered4, strokewidth = 1, color = :orangered4, markercolor = (:white, 0.0), label = L"\text{Hénon map}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 3, 2, 1]; marker = :circle, markersize = 12, strokecolor = :seagreen, strokewidth = 1, color = :seagreen, markercolor = (:white, 0.0), label = L"\text{Lorenz system}")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 4, 2, 1]; marker = :cross, markersize = 12, strokecolor = :violetred4, strokewidth = 1, color = :violetred4, markercolor = (:white, 0.0), label = L"\text{Harmonic oscillator, } x")

hlines!(ax, 0.05; color = :black, linestyle = :dash, label = L"\text{Error}=5%")

axislegend(ax, position = :rt)
# --------------------------------------------------------------------------------
save(plotsdir("series_len_test.png"), fig; px_per_unit = 8)