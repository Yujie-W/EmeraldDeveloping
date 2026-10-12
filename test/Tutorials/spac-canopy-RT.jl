using Test

import Emerald.Namespace as ENS
import Emerald.SPAC as ESPAC
import Emerald.Land as ELAND


@testset "SPAC Level Canopy RT" verbose = true begin
    @testset "Canopy Reflectance Example" begin
        config = ENS.SPACConfig(Float64);
        spac = ENS.BulkSPAC(config);
        ESPAC.initialize_spac!(config, spac);
        ESPAC.spac!(config, spac, 0);
        @test all(spac.canopy.sensor_geometry.auxil.reflectance .> 0);
    end;

    @testset "Change Canopy LAI & SAI" begin
        config = ENS.SPACConfig(Float64);
        spac = ENS.BulkSPAC(config);
        spac.plant.pool.c_pool = Inf;
        ESPAC.prescribe_traits!(config, spac; sai = 0, lai = 3);
        ESPAC.initialize_spac!(config, spac);
        ESPAC.spac!(config, spac, 0);
        @test all(spac.canopy.sensor_geometry.auxil.reflectance .> 0);
        @test ELAND.MODIS_NDVI(config, spac) > 0;
    end;

    @testset "Change Leaf CAB & CAR" begin
        config = ENS.SPACConfig(Float64);
        spac = ENS.BulkSPAC(config);
        spac.plant.pool.c_pool = Inf;
        ESPAC.prescribe_traits!(config, spac; cab = 50.0, car = 6.0);
        ESPAC.initialize_spac!(config, spac);
        ESPAC.spac!(config, spac, 0);
        @test all(spac.canopy.sensor_geometry.auxil.reflectance .> 0);
        @test ELAND.MODIS_NDVI(config, spac) > 0;
    end;

    @testset "Prescribe Soil Albedo" begin
        config = ENS.SPACConfig(Float64);
        config.METHODS.SOIL_ALBEDO = ENS.SoilAlbedoPrescribe()
        spac = ENS.BulkSPAC(config);
        spac.plant.pool.c_pool = Inf;
        spac.soil_bulk.auxil.ρ_sw .= 0.3;
        ESPAC.initialize_spac!(config, spac);
        ESPAC.spac!(config, spac, 3600);
        @test all(spac.canopy.sensor_geometry.auxil.reflectance .> 0);
        @test all(spac.soil_bulk.auxil.ρ_sw .== 0.3);
        @test ELAND.MODIS_NDVI(config, spac) > 0;
    end;
end;
