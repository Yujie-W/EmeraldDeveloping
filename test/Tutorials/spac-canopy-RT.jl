using Test

import Emerald.Namespace as ENS
import Emerald.SPAC as ESPAC
import Emerald.Land as ELAND


@testset "SPAC Level Canopy RT" verbose = true begin
    @testset "Canopy Reflectance" verbose = true begin
        config = ENS.SPACConfig(Float64);
        spac = ENS.BulkSPAC(config);
        ESPAC.initialize_spac!(config, spac);
        spac.plant.pool.c_pool = Inf;
        ESPAC.prescribe_traits!(config, spac; sai = 0, lai = 3);
        ESPAC.spac!(config, spac, 0);
        @test all(spac.canopy.sensor_geometry.auxil.reflectance .> 0);
        @test ELAND.MODIS_NDVI(config, spac) .> 0;
    end;
end;
