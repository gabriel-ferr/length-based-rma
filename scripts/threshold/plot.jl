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
BSON.@load datadir("results_threshold.bson") results
# --------------------------------------------------------------------------------
errors = Array{Float64, 5}(undef, NUM_INITIAL_CONDITIONS, length(THRESHOLD_RANGE), 5, 3, 2)
condensed = Array{Float64, 5}(undef, 3, length(THRESHOLD_RANGE), 5, 3, 2)
for i ∈ eachindex(THRESHOLD_RANGE)
    for j ∈ 1:5
        for k ∈ 1:3
            for n ∈ 1:2
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

mins_vals = Array{Int, 3}(undef, 3, 5, 2)
for n ∈ 1:3
    # --------------------------------------------------------------------------------
    fig = Figure(size = (1200, 1000))
    # --------------------------------------------------------------------------------
    ax = Axis(fig[1, 1], xlabel = L"\text{Threshold}~(\varepsilon)", ylabel = L"\text{Mean relative error}", title = L"\textbf{(A) Gaussian white noise}")
    band!(ax, THRESHOLD_RANGE, condensed[2, :, 1, n, 1], condensed[3, :, 1, n, 1], color = (:blue, 0.2))
    band!(ax, THRESHOLD_RANGE, condensed[2, :, 1, n, 2], condensed[3, :, 1, n, 2], color = (:orange, 0.2))

    lines!(ax, THRESHOLD_RANGE, condensed[1, :, 1, n, 1], label = L"\ell_{min}=1", color = :blue)
    lines!(ax, THRESHOLD_RANGE, condensed[1, :, 1, n, 2], label = L"\ell_{min}=2", color = :orange)

    hlines!(ax, 0.05, linestyle = :dash, color = :black, label = L"\text{Error} = 5\%")
    axislegend(ax, position = :lt)
    # --------------------------------------------------------------------------------
    ax = Axis(fig[1, 2], xlabel = L"\text{Threshold}~(\varepsilon)", ylabel = L"\text{Mean relative error}", title = L"\textbf{(B) Red noise}")
    band!(ax, THRESHOLD_RANGE, condensed[2, :, 2, n, 1], condensed[3, :, 2, n, 1], color = (:blue, 0.2))
    band!(ax, THRESHOLD_RANGE, condensed[2, :, 2, n, 2], condensed[3, :, 2, n, 2], color = (:orange, 0.2))

    lines!(ax, THRESHOLD_RANGE, condensed[1, :, 2, n, 1], label = L"\ell_{min}=1", color = :blue)
    lines!(ax, THRESHOLD_RANGE, condensed[1, :, 2, n, 2], label = L"\ell_{min}=2", color = :orange)

    hlines!(ax, 0.05, linestyle = :dash, color = :black, label = L"\text{Error} = 5\%")
    axislegend(ax, position = :lt)
    # --------------------------------------------------------------------------------
    ax = Axis(fig[2, 1], xlabel = L"\text{Threshold}~(\varepsilon)", ylabel = L"\text{Mean relative error}", title = L"\textbf{(C) Hénon map}")
    band!(ax, THRESHOLD_RANGE, condensed[2, :, 3, n, 1], condensed[3, :, 3, n, 1], color = (:blue, 0.2))
    band!(ax, THRESHOLD_RANGE, condensed[2, :, 3, n, 2], condensed[3, :, 3, n, 2], color = (:orange, 0.2))

    lines!(ax, THRESHOLD_RANGE, condensed[1, :, 3, n, 1], label = L"\ell_{min}=1", color = :blue)
    lines!(ax, THRESHOLD_RANGE, condensed[1, :, 3, n, 2], label = L"\ell_{min}=2", color = :orange)

    hlines!(ax, 0.05, linestyle = :dash, color = :black, label = L"\text{Error} = 5\%")
    axislegend(ax, position = :lt)
    # --------------------------------------------------------------------------------
    ax = Axis(fig[2, 2], xlabel = L"\text{Threshold}~(\varepsilon)", ylabel = L"\text{Mean relative error}", title = L"\textbf{(D) Lorenz system}")
    band!(ax, THRESHOLD_RANGE, condensed[2, :, 4, n, 1], condensed[3, :, 4, n, 1], color = (:blue, 0.2))
    band!(ax, THRESHOLD_RANGE, condensed[2, :, 4, n, 2], condensed[3, :, 4, n, 2], color = (:orange, 0.2))

    lines!(ax, THRESHOLD_RANGE, condensed[1, :, 4, n, 1], label = L"\ell_{min}=1", color = :blue)
    lines!(ax, THRESHOLD_RANGE, condensed[1, :, 4, n, 2], label = L"\ell_{min}=2", color = :orange)

    hlines!(ax, 0.05, linestyle = :dash, color = :black, label = L"\text{Error} = 5\%")
    axislegend(ax, position = :lt)
    # --------------------------------------------------------------------------------
    ax = Axis(fig[3, 1], xlabel = L"\text{Threshold}~(\varepsilon)", ylabel = L"\text{Mean relative error}", title = L"\textbf{(E) Harmonic oscillator, } x")
    band!(ax, THRESHOLD_RANGE, condensed[2, :, 5, n, 1], condensed[3, :, 5, n, 1], color = (:blue, 0.2))
    band!(ax, THRESHOLD_RANGE, condensed[2, :, 5, n, 2], condensed[3, :, 5, n, 2], color = (:orange, 0.2))

    lines!(ax, THRESHOLD_RANGE, condensed[1, :, 5, n, 1], label = L"\ell_{min}=1", color = :blue)
    lines!(ax, THRESHOLD_RANGE, condensed[1, :, 5, n, 2], label = L"\ell_{min}=2", color = :orange)

    hlines!(ax, 0.05, linestyle = :dash, color = :black, label = L"\text{Error} = 5\%")
    axislegend(ax, position = :lt)
    # --------------------------------------------------------------------------------
    save(plotsdir("threshold_n=$n.png"), fig; px_per_unit = 2)
    # --------------------------------------------------------------------------------
    for s ∈ 1:5
        a = [filter(x -> (!isnan(x) && 0 ≤ x ≤ 1), condensed[1, :, s, n, 1]), filter(x -> (!isnan(x)  && 0 ≤ x ≤ 1), condensed[1, :, s, n, 2])]
        mn = [findmin(a[1]), findmin(a[2])]
        mx = [findmax(a[1]), findmax(a[2])]

        println("For n = $n, s = $s : min = $mn; max = $mx")

        mins_vals[n, s, 1] = mn[1][2]
        mins_vals[n, s, 2] = mn[2][2]
    end
    # --------------------------------------------------------------------------------
end
# --------------------------------------------------------------------------------
BSON.@save datadir("min_errors.bson") mins_vals