module LeafLevelSetup

using ..Namespace: SPACCache, SPACConfig, Leaf
using Photosynthesis: C3State, C4State, C3Trait, C4Trait, LeafPhotosystem, LeafPhotosystemAuxil


include("cache.jl");
include("config.jl");
include("leaf.jl");
include("photosystem.jl");


end # module
