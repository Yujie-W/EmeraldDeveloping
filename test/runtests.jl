using Test


@testset verbose = true "Emerald" verbose = true begin
    @testset "ResearchTools" verbose = true begin
        include("ResearchTools/fit-aci.jl");
    end;

    @testset "Tutorials" verbose = true begin
        include("Tutorials/leaf-level-fluorescence.jl");
    end;
end;
