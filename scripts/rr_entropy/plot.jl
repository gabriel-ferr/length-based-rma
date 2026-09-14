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

        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1),


        xgridwidth = 1.0,
        ygridwidth = 1.0,
    )
)
# --------------------------------------------------------------------------------
BSON.@load datadir("results_params.bson") results param_values
# --------------------------------------------------------------------------------
errors = Array{Float64, 5}(undef, NUM_INITIAL_CONDITIONS, length(THRESHOLD_RANGE), 5, 3, 2)
condensed = Array{Float64, 5}(undef, 3, length(THRESHOLD_RANGE), 5, 3, 2)
parameters = Array{Float64, 3}(undef, length(THRESHOLD_RANGE), 5, 2)

for i ∈ eachindex(THRESHOLD_RANGE)
    for j ∈ 1:5
        for n ∈ 1:2
            parameters[i, j, n] = mean(param_values[:, i, j, n])
            for k ∈ 1:3
                for m ∈ 1:NUM_INITIAL_CONDITIONS
                    errors[m, i, j, k, n] = abs(results[m, i, j, k, 2, n] - results[m, i, j, k, 1, n]) / abs(results[m, i, j, k, 1, n])
                    if (isinf(errors[m, i, j, k, n])); errors[m, i, j, k, n] = NaN; end
                    if (errors[m, i, j, k, n] > 1); errors[m, i, j, k, n] = NaN; end
                end

                condensed[1, i, j, k, n] = mean(errors[:, i, j, k, n])
                try
                    condensed[2, i, j, k, n] = quantile(errors[:, i, j, k, n], 0.1)
                    condensed[3, i, j, k, n] = quantile(errors[:, i, j, k, n], 0.9)
                catch
                    condensed[2, i, j, k, n] = NaN
                    condensed[3, i, j, k, n] = NaN
                end
            end
        end
    end
end
# --------------------------------------------------------------------------------
systems = [
    (L"\textbf{(A) Gaussian white noise}", 1, 1), 
    (L"\textbf{(B) Hénon map}", 1, 2),
    (L"\textbf{(C) Lorenz system}", 2, 1),
    (L"\textbf{(D) Harmonic oscillator, } x", 2, 2),
]

pname = [L"\text{Recurrence rate}~(RR)", L"\text{Recurrence microstates entropy}~(RME)"]
# --------------------------------------------------------------------------------
for k ∈ 1:3
    for n ∈ 1:2
        fig = Figure(size = (1200, 680))

        for s ∈ eachindex(systems)
            ax = Axis(fig[systems[s][2], systems[s][3]], xlabel = pname[n], ylabel = L"\text{Relative error}", title = systems[s][1], yscale = log10)
            idx = sortperm(parameters[:, s, n])

            if (n == 1)
                vspan!(ax, 0.1, 0.8, color = (:gray, 0.15))
            end

            x = parameters[idx, s, n]
            y = condensed[1, idx, s, k, :]
            error = condensed[2:3, idx, s, k, :]
            band!(ax, x, error[1, :, 1], error[2, :, 1], color = (:blue, 0.08))
            band!(ax, x, error[1, :, 2], error[2, :, 2], color = (:orangered4, 0.08))

            scatterlines!(ax, x, y[:, 1], label = L"\ell_{min}=1", color = :blue, marker = :utriangle, markersize = 8, strokecolor = :blue, strokewidth = 1, markercolor = (:white, 0.0), linewidth = 0.8)
            scatterlines!(ax, x, y[:, 2], label = L"\ell_{min}=2", color = :orangered4, marker = :rect, markersize = 8, strokecolor = :orangered4, strokewidth = 1, markercolor = (:white, 0.0), linewidth = 0.8)

            hlines!(ax, 0.05, linestyle = :dash, color = :black, label = L"\text{Error} = 5\%")

            axislegend(ax, position = n == 1 ? :rt : :lt)
        end
        save(plotsdir("params_p=$(n)_k=$(k).png"), fig; px_per_unit = 6)
        save(plotsdir("params_p=$(n)_k=$(k).pdf"), fig)
    end
end
# --------------------------------------------------------------------------------