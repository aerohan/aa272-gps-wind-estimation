module GPSModel

using Dates
using PyCall
using NLsolve
using LinearAlgebra
using StaticArrays
using Geodesy

const DATADIR = joinpath(@__DIR__, "..", "..", "data")

include("ephemeris.jl")
include("meas.jl")

export  broadcast_ephemeris,
        EphemerisErrors,
        compute_pseudoranges_truth!,
        compute_pseudoranges_nav!

end # module
