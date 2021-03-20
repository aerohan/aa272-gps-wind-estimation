module UAVWindEstimation

using DynamicsAndControl
using StaticArrays
using LinearAlgebra
using Dates
using Geodesy
using GPSModel
using Plots
using GMT: coast, arrows!
using DataFrames
using Interpolations
using CSV
using Measures: mm

const DATADIR = joinpath(@__DIR__, "..", "..", "data")

include("models/atmo.jl")
include("estimation/ukf.jl")
include("sim/static_gnss.jl")
include("sim/uav_flight.jl")

export  simulate,   # from DynamicsAndControl
        plot        # from RecipesBase (Plots)

end # module
