# function canopy_turbulence!(dt, physcon, forcvar, surfvar, leafvar, soilvar, fluxvar)




function canopy_turbulence!(spac::BulkSPAC{FT}, dt::Number) where {FT<:AbstractFloat}
    # define the constants used for the calculations
    k_von_karman = K_VON_KARMAN(FT);
    drag_coef = FT(0.25);

    # Canopy density length scale
    Lc::FT = spac.plant.zs[end] / (drag_coef * (spac.canopy.lai + spac.canopy.sai));

    #=
	obu0 = 100.0
	obu1 = -100.0
	tol = 0.01
	fluxvar, _ = hybrid_root("obukhov_function", physcon, forcvar, surfvar, fluxvar, obu0, obu1, tol)
	fluxvar.obu[p] = fluxvar.obu_ustar[p]

	h_minus_d = surfvar.hc[p] - fluxvar.disp[p]
	psi_m_hc = psi_m_monin_obukhov(h_minus_d / fluxvar.obu[p])
	psi_m_rsl_hc = psi_m_rsl(h_minus_d, h_minus_d, fluxvar.obu[p], fluxvar.c1m[p], fluxvar.c2)

	for ic in (surfvar.ntop[p] + 1):surfvar.nlev[p]
		z_minus_d = surfvar.zs[p, ic] - fluxvar.disp[p]
		psi_m_zs = psi_m_monin_obukhov(z_minus_d / fluxvar.obu[p])
		psi_m_rsl_zs = psi_m_rsl(z_minus_d, h_minus_d, fluxvar.obu[p], fluxvar.c1m[p], fluxvar.c2)
		psim = -psi_m_zs + psi_m_hc + psi_m_rsl_zs - psi_m_rsl_hc + k_von_karman / fluxvar.beta[p]
		fluxvar.wind[p, ic] = fluxvar.ustar[p] / k_von_karman * (log(z_minus_d / h_minus_d) + psim)
	end

	fluxvar.uaf[p] = fluxvar.ustar[p] / fluxvar.beta[p]
	lm = 2 * fluxvar.beta[p]^3 * fluxvar.Lc[p]
	lm_over_beta = lm / fluxvar.beta[p]
	for ic in (surfvar.nsoi[p] + 1):surfvar.ntop[p]
		fluxvar.wind[p, ic] = fluxvar.uaf[p] * exp((surfvar.zs[p, ic] - surfvar.hc[p]) / lm_over_beta)
		fluxvar.wind[p, ic] = max(fluxvar.wind[p, ic], 0.1)
	end
	fluxvar.wind[p, surfvar.nsoi[p]] = 0

	h_minus_d = surfvar.hc[p] - fluxvar.disp[p]
	psi_c_hc = psi_c_monin_obukhov(h_minus_d / fluxvar.obu[p])
	psi_c_rsl_hc = psi_c_rsl(h_minus_d, h_minus_d, fluxvar.obu[p], fluxvar.c1c[p], fluxvar.c2)

	for ic in (surfvar.ntop[p] + 1):(surfvar.nlev[p] - 1)
		z_minus_d = surfvar.zs[p, ic] - fluxvar.disp[p]
		psi_c_zs = psi_c_monin_obukhov(z_minus_d / fluxvar.obu[p])
		psi_c_rsl_zs = psi_c_rsl(z_minus_d, h_minus_d, fluxvar.obu[p], fluxvar.c1c[p], fluxvar.c2)
		psic1 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

		z_minus_d = surfvar.zs[p, ic + 1] - fluxvar.disp[p]
		psi_c_zs = psi_c_monin_obukhov(z_minus_d / fluxvar.obu[p])
		psi_c_rsl_zs = psi_c_rsl(z_minus_d, h_minus_d, fluxvar.obu[p], fluxvar.c1c[p], fluxvar.c2)
		psic2 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

		psic = psic2 - psic1
		zlog = log((surfvar.zs[p, ic + 1] - fluxvar.disp[p]) / (surfvar.zs[p, ic] - fluxvar.disp[p]))
		fluxvar.ga_prof[p, ic] = forcvar.rhomol[p] * k_von_karman * fluxvar.ustar[p] / (zlog + psic)
	end

	ic = surfvar.nlev[p]
	z_minus_d = surfvar.zs[p, ic] - fluxvar.disp[p]
	psi_c_zs = psi_c_monin_obukhov(z_minus_d / fluxvar.obu[p])
	psi_c_rsl_zs = psi_c_rsl(z_minus_d, h_minus_d, fluxvar.obu[p], fluxvar.c1c[p], fluxvar.c2)
	psic1 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

	z_minus_d = forcvar.zref[p] - fluxvar.disp[p]
	psi_c_zs = psi_c_monin_obukhov(z_minus_d / fluxvar.obu[p])
	psi_c_rsl_zs = psi_c_rsl(z_minus_d, h_minus_d, fluxvar.obu[p], fluxvar.c1c[p], fluxvar.c2)
	psic2 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

	psic = psic2 - psic1
	zlog = log((forcvar.zref[p] - fluxvar.disp[p]) / (surfvar.zs[p, ic] - fluxvar.disp[p]))
	fluxvar.ga_prof[p, ic] = forcvar.rhomol[p] * k_von_karman * fluxvar.ustar[p] / (zlog + psic)

	for ic in (surfvar.nsoi[p] + 1):(surfvar.ntop[p] - 1)
		zl = surfvar.zs[p, ic] - surfvar.hc[p]
		zu = surfvar.zs[p, ic + 1] - surfvar.hc[p]
		res = fluxvar.PrSc[p] / (fluxvar.beta[p] * fluxvar.ustar[p]) *
			  (exp(-zl / lm_over_beta) - exp(-zu / lm_over_beta))
		fluxvar.ga_prof[p, ic] = forcvar.rhomol[p] / res
	end

	ic = surfvar.ntop[p]
	zl = surfvar.zs[p, ic] - surfvar.hc[p]
	zu = 0
	res = fluxvar.PrSc[p] / (fluxvar.beta[p] * fluxvar.ustar[p]) *
		  (exp(-zl / lm_over_beta) - exp(-zu / lm_over_beta))
	ga_below_hc = forcvar.rhomol[p] / res

	z_minus_d = surfvar.hc[p] - fluxvar.disp[p]
	psi_c_zs = psi_c_monin_obukhov(z_minus_d / fluxvar.obu[p])
	psi_c_rsl_zs = psi_c_rsl(z_minus_d, h_minus_d, fluxvar.obu[p], fluxvar.c1c[p], fluxvar.c2)
	psic1 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

	z_minus_d = surfvar.zs[p, ic + 1] - fluxvar.disp[p]
	psi_c_zs = psi_c_monin_obukhov(z_minus_d / fluxvar.obu[p])
	psi_c_rsl_zs = psi_c_rsl(z_minus_d, h_minus_d, fluxvar.obu[p], fluxvar.c1c[p], fluxvar.c2)
	psic2 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

	psic = psic2 - psic1
	zlog = log((surfvar.zs[p, ic + 1] - fluxvar.disp[p]) / (surfvar.hc[p] - fluxvar.disp[p]))
	ga_above_hc = forcvar.rhomol[p] * k_von_karman * fluxvar.ustar[p] / (zlog + psic)
	fluxvar.ga_prof[p, ic] = 1 / (1 / ga_below_hc + 1 / ga_above_hc)

	sumres = 1 / ga_above_hc
	for ic in (surfvar.ntop[p] + 1):surfvar.nlev[p]
		sumres += 1 / fluxvar.ga_prof[p, ic]
	end
	if abs(1 / sumres - fluxvar.gac[p]) > 1e-6
		error("canopy_turbulence: above-canopy aerodynamic conductance error")
	end

	ic = surfvar.nsoi[p]
	ic_plus_one = ic + 1
	z0m_g = 0.01
	z0c_g = 0.1 * z0m_g
	zlog_m = log(surfvar.zs[p, ic_plus_one] / z0m_g)
	zlog_c = log(surfvar.zs[p, ic_plus_one] / z0c_g)
	ustar_g = fluxvar.wind[p, ic_plus_one] * k_von_karman / zlog_m
	ustar_g = max(ustar_g, 0.01)
	res = zlog_c / (k_von_karman * ustar_g)
	fluxvar.ga_prof[p, ic] = forcvar.rhomol[p] / res

	for ic in surfvar.nsoi[p]:surfvar.nlev[p]
		res = min(forcvar.rhomol[p] / fluxvar.ga_prof[p, ic], 500)
		fluxvar.ga_prof[p, ic] = forcvar.rhomol[p] / res
	end

	fluxvar = scalar_profile(dt, p, physcon, forcvar, surfvar, leafvar, soilvar, fluxvar)
	fluxvar.taf[p] = fluxvar.tair[p, surfvar.ntop[p]]
	fluxvar.qaf[p] = fluxvar.qair[p, surfvar.ntop[p]]
	return fluxvar
    =#
end
