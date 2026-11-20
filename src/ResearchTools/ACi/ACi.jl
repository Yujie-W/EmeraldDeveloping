module ACi

using DataFrames: DataFrame
using PkgUtility.MathTools: ReduceStepMethodND, SolutionToleranceND, find_peak, nanmax, nanmin, rmse

using Photosynthesis: LeafPhotosystem
using Photosynthesis: C3Trait, C4Trait
using Photosynthesis: AcMethodC3Vcmax, AcMethodC4Vcmax, AjMethodC3Jmax, AjMethodC3Vqmax, AjMethodC4JPSII, ApMethodC3Vcmax, ApMethodC4Vcmax, ApMethodC4Vpmax
using Photosynthesis: photosynthesis!, temperature_dependence!, temperature_correction

using ..Namespace: SPACConfig
using ..Namespace: AirLayer

using ..LeafLevelSetup: leaf_level_photosystem, leaf_level_spac_cache


include("curve.jl");
include("rmse.jl");
include("fit.jl");

include("procedure/c3-vcmax-jmax-vcmax.jl");
include("procedure/c3-vcmax-vqmax-vcmax.jl");
include("procedure/c4-vcmax-jpsii-vcmaxpi.jl");
include("procedure/c4-vcmax-jpsii-vpmaxpi.jl");

include("pipeline.jl");


end # module
