module CanopyAirSpace

using PkgUtility.MathTools: NewtonBisectionMethod, SolutionTolerance, find_zero, upper_quadratic
using PkgUtility.UniversalConstants: CP_D, GRAVITY, K_VON_KARMAN, R_V

using ..Namespace: BulkSPAC


include("canopy_turbulence.jl");
include("flux–gradient.jl");
include("refactor.jl");
include("xiaomiao.jl");
include("wind-ratio.jl");


end # module CanopyAirSpace
