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
BSON.@load datadir("results_ho.bson") results
# --------------------------------------------------------------------------------
errors = Array{Float64, 3}(undef, NUM_INITIAL_CONDITIONS, length(TIMESERIES_LENGTH_RANGE), 3)
condensed = Array{Float64, 3}(undef, length(TIMESERIES_LENGTH_RANGE), 3, 3)

for i ∈ eachindex(TIMESERIES_LENGTH_RANGE)
    for j ∈ 1:3
        for k ∈ 1:NUM_INITIAL_CONDITIONS
            errors[k, i, j] = abs(results[k, i, j, 2] - results[k, i, j, 1]) / abs(results[k, i, j, 1])
            if (isinf(errors[k, i, j])); errors[k, i, j] = NaN; end
            if (errors[k, i, j] > 1); errors[k, i, j] = NaN; end
        end

        condensed[i, j, 1] = mean(errors[:, i, j])
        try
            condensed[i, j, 2] = quantile(errors[:, i, j], 0.1)
            condensed[i, j, 3] = quantile(errors[:, i, j], 0.9)
        catch
            condensed[i, j, 2] = NaN
            condensed[i, j, 3] = NaN
        end
    end
end
# --------------------------------------------------------------------------------
fig = Figure(size = (620, 465))

ax = Axis(fig[1, 1], xlabel = L"\text{Time series length}", ylabel = L"\text{Relative error}")
band!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 1, 2], condensed[:, 1, 3]; color = (:blue, 0.1))
band!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 2, 2], condensed[:, 2, 3]; color = (:orangered4, 0.1))

scatterlines!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 1, 1]; marker = :utriangle, markersize = 12, strokecolor = :blue, strokewidth = 1, color = :blue, markercolor = (:white, 0.0), label = L"\langle D_2 \rangle")
scatterlines!(ax, TIMESERIES_LENGTH_RANGE, condensed[:, 2, 1]; marker = :rect, markersize = 12, strokecolor = :orangered4, strokewidth = 1, color = :orangered4, markercolor = (:white, 0.0), label = L"\langle TT_2 \rangle")

hlines!(ax, 0.05; color = :black, linestyle = :dash, label = L"\text{Error}=5%")

axislegend(ax, position = :rt)
# --------------------------------------------------------------------------------
save(plotsdir("harmonic_test.png"), fig; px_per_unit = 4)