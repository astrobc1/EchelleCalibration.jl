module EchelleCalibration

using Pkg

using NaNStatistics
using LsqFit
using Infiltrator
using FindPeaks1D
using DataInterpolations
using Polynomials
using ImageFiltering
using PyCall
using StatsBase

using Pkg
Pkg.develop(path="/Users/cale/Codes/JuliaProjects/Echelle/")

const SPEED_OF_LIGHT_MPS = 299792458.0


include("utils.jl")
include("doppler.jl")
include("mode_finding.jl")
include("mode_fitting.jl")
include("pixels.jl")
include("lfc.jl")
include("poly2d.jl")
include("drifts.jl")
include("lsf.jl")

end
