# --------------------------------------------------------------------------------
#       Global settings file
# --------------------------------------------------------------------------------

#   [RECURRENCE SETTINGS]
const THRESHOLD_RANGE = range(0.0, 1.0, 100)

const THRESHOLD_RANGE_SAMPLING = range(0.01, 0.5, 25)
const SAMPLING_RATIO_RANGE = [0.01, 0.05, 0.1, 0.25, 0.5]

const FIXED_THRESHOLD = 0.5

#   [TIME SERIES SETTINGS]
const NUM_INITIAL_CONDITIONS = 100

const TRANSIENT_TIME = 10000
const DEFAULT_TIMESERIES_LENGTH = 3000
const TIMESERIES_LENGTH_RANGE = ceil.(Int, collect(range(1000, 10000, 40)))

#   Hénon
const henon_params = [1.4, 0.3]
#   Lorenz
const lorenz_params = [10.0, 8/3, 28.0]
#   Harmonic oscillator
const ho_params = 2π

# --------------------------------------------------------------------------------
#       DON'T CHANGE!
# --------------------------------------------------------------------------------
mkpath(datadir())
mkpath(datadir("raw"))

mkpath(plotsdir())
# --------------------------------------------------------------------------------