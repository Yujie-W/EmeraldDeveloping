#
# These tests requires downloading a lot of data, and thus we only test them locally
#
using Test


@testset verbose = true "Emerald" verbose = true begin
    @testset "Continous Integration" verbose = true begin
        include("Land/site-run.jl");
    end;
end;
