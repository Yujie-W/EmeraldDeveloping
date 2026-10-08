using PkgUtility.DataIO: read_csv
using Test

import Emerald.Namespace as ENS
import Emerald.ResearchTools as ERT
import Photosynthesis as PS


@testset "Emerald ResearchTools" verbose = true begin
    df3 = read_csv(joinpath(@__DIR__, "../..", "data/examples", "C3-ACi.csv"));
    df4 = read_csv(joinpath(@__DIR__, "../..", "data/examples", "C4-ACi.csv"));
    df3.T_LEAF .+= 273.15;  # convert to Kelvin
    df4.T_LEAF .+= 273.15;  # convert to Kelvin

    @testset "C3 Jmax" begin
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.PS_METHODS.C3_AC_METHOD = PS.AcMethodC3VcmaxPi();
        config.METHODS.PS_METHODS.C3_AJ_METHOD = PS.AjMethodC3JmaxPi();
        config.METHODS.PS_METHODS.C3_AP_METHOD = PS.ApMethodC3Vcmax();
        config.METHODS.PS_METHODS.COLIMIT_J = PS.ColimitJCLM(Float64);
        config.METHODS.PS_METHODS.FLUORESCENCE_METHOD_C3 = PS.KNFluorescenceModel{Float64}();
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "Jmax25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "Jmax25", "Rd25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "Jmax25", "Γstar25", "Rd25"]);
        @test !any(isnan.(result[1]));
    end;

    @testset "C3 Vqmax" begin
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.PS_METHODS.C3_AC_METHOD = PS.AcMethodC3VcmaxPi();
        config.METHODS.PS_METHODS.C3_AJ_METHOD = PS.AjMethodC3VqmaxPi();
        config.METHODS.PS_METHODS.C3_AP_METHOD = PS.ApMethodC3Vcmax();
        config.METHODS.PS_METHODS.COLIMIT_J = PS.SerialColimit();
        config.METHODS.PS_METHODS.FLUORESCENCE_METHOD_C3 = PS.CytochromeFluorescenceModel();
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "b₆f"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "b₆f", "Rd25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "b₆f", "Γstar25", "Rd25"]);
        @test !any(isnan.(result[1]));
    end;

    @testset "C4 Vcmax" begin
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.PS_METHODS.C4_AP_METHOD = PS.ApMethodC4VcmaxPi();
        config.METHODS.PS_METHODS.FLUORESCENCE_METHOD_C4 = PS.KNFluorescenceModel{Float64}();
        result = ERT.ACi.aci_fit!(config, df4, "C4", ["Vcmax25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df4, "C4", ["Vcmax25", "Rd25"]);
        @test !any(isnan.(result[1]));
    end;

    @testset "C4 Vpmax" begin
        config = ERT.LeafLevelSetup.leaf_level_config(Float64);
        config.METHODS.PS_METHODS.C4_AP_METHOD = PS.ApMethodC4VpmaxPi();
        config.METHODS.PS_METHODS.FLUORESCENCE_METHOD_C4 = PS.KNFluorescenceModel{Float64}();
        result = ERT.ACi.aci_fit!(config, df4, "C4", ["Vcmax25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df4, "C4", ["Vcmax25", "Vpmax25"]);
        @test !any(isnan.(result[1]));
        result = ERT.ACi.aci_fit!(config, df4, "C4", ["Vcmax25", "Vpmax25", "Rd25"]);
        @test !any(isnan.(result[1]));
    end;
end;
