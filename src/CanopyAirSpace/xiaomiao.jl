# using QuadGK: quadgk

const RSL_C2 = 0.5

const DEFAULT_UREF = 5.169;
const DEFAULT_MIN_WIND = 0.1;
const RD_AIR = 287.05;
const P0 = 100000.0;

function checked_number(name, value; positive=false, nonnegative=false, finite=true)
    value isa Real || throw(ArgumentError("$name must be a real scalar."))
    x = Float64(value)
    isnan(x) && throw(ArgumentError("$name must not be NaN."))
    finite && !isfinite(x) && throw(ArgumentError("$name must be finite."))
    positive && x <= 0 && throw(ArgumentError("$name must be positive."))
    nonnegative && x < 0 && throw(ArgumentError("$name must be nonnegative."))
    return x
end




# Integrate the COMPLETE gradient over a finite interval. This is algebraically
# identical to log - psi + psi_hat differences in GetPsiRSL, without subtracting
# two nearly equal improper integrals. t=(z-d)/(hc-d) is dimensionless.
function profile_factor(z, state; scalar=false)
    z >= state.hc || throw(ArgumentError("RSL factors require z >= canopy height."))
    offset = scalar ? 0.0 : K_VON_KARMAN()/state.β
    z == state.hc && return offset
    tmax = max(1.0,(z-state.d)/state.hd)
    tmax == 1 && return offset
    c1 = scalar ? state.c₁c : state.c₁m
    phi = scalar ? Φ_heat_vapor : Φ_momentum
    integrand(t) = phi(t*state.hd*state.invL) *
                   (1-c1*exp(-RSL_C2*t/2))/t
    val, _ = quadgk(integrand, 1.0, tmax; rtol=1e-9, atol=1e-11)
    factor = offset + val
    isfinite(factor) && factor > 0 ||
        error("Invalid RSL factor; check geometry and stability.")
    return factor
end


function measured_wind(uref, wind_speed, wind_u, wind_v, wind_height, zref;
                       allow_demo_default=true)
    components = wind_u !== nothing || wind_v !== nothing
    ninputs = (uref !== nothing) + (wind_speed !== nothing) + components
    ninputs <= 1 || throw(ArgumentError("Supply only uref, wind_speed, or wind_u/wind_v."))
    if uref !== nothing
        wind_height === nothing || wind_height == zref ||
            throw(ArgumentError("uref is already at zref; use wind_speed for another height."))
        return checked_number("uref", uref; nonnegative=true), zref, :uref
    elseif components
        wind_u !== nothing && wind_v !== nothing ||
            throw(ArgumentError("Both wind_u and wind_v are required (horizontal m/s)."))
        speed = hypot(checked_number("wind_u", wind_u), checked_number("wind_v", wind_v))
        source = :components
    elseif wind_speed !== nothing
        speed = checked_number("wind_speed", wind_speed; nonnegative=true)
        source = :speed
    else
        allow_demo_default || throw(ArgumentError("Supply measured uref, wind_speed, or wind_u/wind_v; SPAC wind fields are not inferred."))
        wind_height === nothing || throw(ArgumentError("wind_height requires a measured wind input."))
        return DEFAULT_UREF, zref, :demo_default
    end
    zwind = wind_height === nothing ? zref :
            checked_number("wind_height", wind_height; positive=true)
    return speed, zwind, source
end

# Solve the source's OUTER zeta bound with geometry recomputed at every trial.
# The near-neutral |zeta|>=0.01 clamp is deliberately omitted: L=Inf is exact.
function inverse_limit(signL, hc, pai, zbound)
    target = signL > 0 ? 1.0 : 2.0
    lo, hi = 0.0, 1/(zbound-hc)
    while (zbound-canopy_state(hc,pai,signL*hi).d)*hi < target
        hi *= 2
    end
    for _ in 1:70
        mid = (lo+hi)/2
        if (zbound-canopy_state(hc,pai,signL*mid).d)*mid > target
            hi = mid
        else
            lo = mid
        end
    end
    return (lo+hi)/2
end

function prescribed_state(hc, pai, zbound, L, stability_limits)
    L = checked_number("L", L; finite=false)
    L != 0 || throw(ArgumentError("L must be nonzero; Inf denotes neutral conditions."))
    invL = 1/L
    limited = false
    if invL != 0 && stability_limits == :source
        limit = inverse_limit(sign(invL),hc,pai,zbound)
        limited = abs(invL) > limit
        invL = sign(invL)*min(abs(invL),limit)
    end
    return canopy_state(hc,pai,invL), L, limited
end

function required_fields(data, names, label)
    data isa NamedTuple || throw(ArgumentError("$label must be a NamedTuple."))
    missing_names = filter(n -> !hasproperty(data,n), names)
    isempty(missing_names) || throw(ArgumentError("$label is missing: $(join(missing_names, ", "))."))
end

function flux_state(; tair, pair, qair, H, LE, cp=CP_D(), latent_heat=nothing)
    tair = checked_number("tair (K)",tair; positive=true)
    pair = checked_number("pair (Pa)",pair; positive=true)
    qair = checked_number("qair (kg/kg)",qair; nonnegative=true)
    qair < 1 || throw(ArgumentError("qair must be < 1 kg/kg, not RH or g/kg."))
    H = checked_number("H (W/m^2)",H)
    LE = checked_number("LE (W/m^2)",LE)
    cp = checked_number("cp (J/kg/K)",cp; positive=true)
    Lv = latent_heat === nothing ? 2.501e6-2361*(tair-273.15) :
         checked_number("latent_heat (J/kg)",latent_heat; positive=true)
    Lv > 0 || throw(ArgumentError("Temperature gives nonpositive latent heat."))
    rho = pair/(RD_AIR*tair*(1+(R_V() /RD_AIR-1)*qair))
    exner_inverse = (P0/pair)^(RD_AIR/cp)
    theta_ref = tair*exner_inverse
    theta_v = theta_ref*(1+0.61*qair)
    Ftheta = H/(rho*cp)*exner_inverse
    Fq = LE/(rho*Lv)
    # Same linearized virtual-potential-temperature scale as source ObuFunc.
    Fvirtual = Ftheta + 0.61*theta_ref*Fq
    return (; tair,pair,qair,H,LE,cp,Lv,rho,theta_ref,theta_v,Ftheta,Fq,Fvirtual)
end

"""
    obukhov_length(; ustar, tair, pair, qair, H, LE,
                     cp=1004.64, latent_heat=nothing)

Compute physical L [m] from known friction velocity [m/s] and ABOVE-CANOPY
fluxes. tair [K], pair [Pa], qair [kg/kg]; H and LE [W/m^2], positive upward.
The source's approximation is Fv = Ftheta + 0.61*theta_ref*Fq.
Neutral Fv=0 gives Inf. Nonzero Fv with ustar=0 is outside this shear model.
This helper does not impose the optional profile stability bounds.
"""
function obukhov_length(; ustar, kwargs...)
    ustar = checked_number("ustar",ustar; nonnegative=true)
    therm = flux_state(; kwargs...)
    therm.Fvirtual == 0 && return Inf
    ustar > 0 || throw(ArgumentError("Nonzero buoyancy flux requires nonzero ustar in this model."))
    return -therm.theta_v*ustar^3/(K_VON_KARMAN()*GRAVITY()*therm.Fvirtual)
end

function gradient_state(gradients)
    required_fields(gradients,(:theta_ref,:theta_canopy,:q_ref,:q_canopy),"gradients")
    tr = checked_number("theta_ref (K)",gradients.theta_ref; positive=true)
    tc = checked_number("theta_canopy (K)",gradients.theta_canopy; positive=true)
    qr = checked_number("q_ref (kg/kg)",gradients.q_ref; nonnegative=true)
    qc = checked_number("q_canopy (kg/kg)",gradients.q_canopy; nonnegative=true)
    qr < 1 && qc < 1 || throw(ArgumentError("Specific humidities must be < 1 kg/kg."))
    return (; theta_ref=tr, theta_canopy=tc, q_ref=qr, q_canopy=qc,
              theta_v=tr*(1+0.61*qr), delta_virtual=tr-tc+0.61*tr*(qr-qc))
end

# Near critical stable fluxes, two roots can be much closer than a scan step.
# Refine each detected local maximum before declaring that no root exists.
function residual_maximum(f, left, right)
    ratio = (sqrt(5.0)-1)/2
    x1 = right-ratio*(right-left)
    x2 = left+ratio*(right-left)
    f1,f2 = f(x1),f(x2)
    for _ in 1:80
        if f1 < f2
            left,x1,f1 = x1,x2,f2
            x2 = left+ratio*(right-left)
            f2 = f(x2)
        else
            right,x2,f2 = x2,x1,f1
            x1 = right-ratio*(right-left)
            f1 = f(x1)
        end
    end
    return f1 >= f2 ? (x1,f1) : (x2,f2)
end

function automatic_state(hc,pai,zref,speed,zwind,zbound,fluxes,gradients,
                         stability_limits,solver_rtol,maxiter)
    (fluxes === nothing) != (gradients === nothing) ||
        throw(ArgumentError("L=:auto requires exactly one of fluxes or gradients."))
    if fluxes !== nothing
        required_fields(fluxes,(:tair,:pair,:qair,:H,:LE),"fluxes")
        therm = flux_state(; fluxes...)
        source = :fluxes
        driving = -therm.Fvirtual
    else
        therm = gradient_state(gradients)
        source = :gradients
        driving = therm.delta_virtual
    end
    function evaluate(invL)
        s = canopy_state(hc,pai,invL)
        ustar = speed*K_VON_KARMAN()/profile_factor(zwind,s)
        if source == :fluxes
            target = -K_VON_KARMAN()*GRAVITY()*therm.Fvirtual/(therm.theta_v*ustar^3)
        else
            Fc = profile_factor(zref,s; scalar=true)
            tvstar = K_VON_KARMAN()*therm.delta_virtual/Fc
            target = K_VON_KARMAN()*GRAVITY()*tvstar/(therm.theta_v*ustar^2)
        end
        return s,target
    end
    if driving == 0
        return canopy_state(hc,pai,0.0),
               (; mode=source,status=:neutral,closure_converged=true,iterations=0,
                  L_physical=Inf,limited=false,inverse_L_residual=0.0,thermal=therm)
    end
    speed > 0 || throw(ArgumentError("Nonzero buoyancy forcing with zero measured wind needs a free-convection model."))
    signL = sign(driving)
    # Without source bounds, search down to |L|=0.1 m, then fail explicitly.
    cap = stability_limits == :source ? inverse_limit(signL,hc,pai,zbound) : 10.0
    _, target0 = evaluate(0.0)
    start = min(abs(target0)/1024,cap/1e7)
    start = max(start,nextfloat(0.0))
    lo, flo, hi = 0.0, -abs(target0), start
    previous, fprevious = lo,flo
    iterations = 0
    bracketed = false
    while true
        _, target = evaluate(signL*hi)
        fhi = hi-abs(target)
        if fhi >= 0
            bracketed = true
            break
        end
        if fhi < flo && flo > fprevious
            residual(x) = x-abs(last(evaluate(signL*x)))
            peak,fpeak = residual_maximum(residual,previous,hi)
            if fpeak >= 0
                lo,hi = previous,peak
                bracketed = true
                break
            end
        end
        previous,fprevious = lo,flo
        lo, flo = hi, fhi
        hi >= cap && break
        hi = min(max(hi*1.25,nextfloat(hi)),cap)  # First (weak) branch from neutrality.
    end
    if bracketed
        for i in 1:maxiter
            iterations = i
            mid = (lo+hi)/2
            _,target = evaluate(signL*mid)
            fmid = mid-abs(target)
            if abs(fmid) <= solver_rtol*max(mid,abs(target),floatmin(Float64))
                lo = hi = mid
                break
            end
            if fmid >= 0
                hi = mid
            else
                lo, flo = mid, fmid
            end
        end
        invL = signL*(lo+hi)/2
    elseif stability_limits == :source
        invL = signL*cap
    else
        error("No self-consistent L found for |L| >= 0.1 m. Stable prescribed-flux cases can have no weak-branch solution; check inputs or use stability_limits=:source for an explicitly flagged bound.")
    end
    s,target = evaluate(invL)
    residual = invL-target
    converged = abs(residual) <= 2*solver_rtol*max(abs(invL),abs(target),floatmin(Float64))
    bracketed && !converged && error("Obukhov solver failed to converge within maxiter=$maxiter.")
    diagnostic = (; mode=source,status=bracketed ? :converged : :limited,
                   closure_converged=converged,iterations,L_physical=1/target,
                   limited=!bracketed,inverse_L_residual=residual,thermal=therm)
    return s,diagnostic
end

"""
    reference_wind_speed(; wind_speed=nothing, wind_u=nothing, wind_v=nothing,
                           wind_height, zref, hc, pai, L=Inf,
                           stability_limits=:source)

Horizontal wind at zref [m/s]. Components are combined with hypot(u,v).
For wind_height != zref use the SAME canopy/RSL model, not a power law.
Both heights must be above the canopy. L must be prescribed here; use
wind_profile(..., L=:auto) to couple the height conversion to thermal forcing.
"""
function reference_wind_speed(; wind_speed=nothing,wind_u=nothing,wind_v=nothing,
                                wind_height,zref,hc,pai,L=Inf,stability_limits=:source)
    stability_limits in (:source,:none) || throw(ArgumentError("stability_limits must be :source or :none."))
    speed,zwind,_ = measured_wind(nothing,wind_speed,wind_u,wind_v,wind_height,zref;
                                 allow_demo_default=false)
    zwind > hc || throw(ArgumentError("wind_height must exceed hc."))
    s,_,_ = prescribed_state(hc,pai,max(zref,zwind),L,stability_limits)
    return speed*profile_factor(zref,s)/profile_factor(zwind,s)
end

function height_grid(hc,zref,dz,heights)
    if heights !== nothing
        heights isa AbstractVector || throw(ArgumentError("heights must be a vector in m."))
        zz = [checked_number("height",x; nonnegative=true) for x in heights]
        all(x -> x <= zref,zz) || throw(ArgumentError("heights must not exceed zref."))
        return sort!(unique(vcat(0.0,zz,hc,zref)))
    end
    nbelow = max(1,round(Int,hc/dz,RoundNearestTiesAway))
    nabove = max(1,round(Int,(zref-hc)/dz,RoundNearestTiesAway))
    below = collect(range(0.0,hc; length=nbelow+1))
    above = collect(range(hc,zref; length=nabove+1))
    return sort!(unique(vcat(0.0,(below[1:end-1]+below[2:end])/2,hc,
                            (above[1:end-1]+above[2:end])/2,zref)))
end

"""
    wind_profile(; hc=21, pai=5.051612734794617, zref=46, uref=nothing,
                   wind_speed=nothing, wind_u=nothing, wind_v=nothing,
                   wind_height=nothing, dz=0.5, heights=nothing, L=Inf,
                   fluxes=nothing, gradients=nothing, minWind=0.1,
                   stability_limits=:source, solver_rtol=1e-8, maxiter=100,
                   outputDir=joinpath(@__DIR__, "wind_profile_results"),
                   saveResults=true, makePlot=true, showPlot=true, verbose=true)

uref is already at zref; alternatively give measured wind_speed or both wind
components, with wind_height (defaults to zref). No-input scalar calls retain
the original DEMO wind 5.169 m/s. wind_profile(spac; ...) requires real wind input.

L=Inf: prescribed neutral; L::Real: prescribed stability; L=:auto: solve L and
ustar together using exactly one thermal NamedTuple:
  fluxes=(;tair, pair, qair, H, LE)  # K, Pa, kg/kg, W/m^2 upward positive
  gradients=(;theta_ref, theta_canopy, q_ref, q_canopy) # K, K, kg/kg, kg/kg
Both gradient temperatures are POTENTIAL temperatures, at zref and hc.
Fluxes must describe the above-canopy turbulent exchange, not a single leaf.
Thermal inputs with prescribed L are rejected to prevent accidental omission.

stability_limits=:source enforces -2 <= (max(zref,wind_height)-d)/L <= 1,
while retaining exact neutrality; :none leaves L unbounded in prescribed mode.
A bound can prevent thermal closure: inspect diagnostics.closure_converged,
diagnostics.L_physical, config.L, and diagnostics.inverse_L_residual.
This is the supplied Fortran/v1 parameterization, not the 2026 v2 revision.

Wind_unfloored_m_s is the model wind. Wind_m_s applies a canopy transport floor
min(minWind,U(hc)) at positive subcanopy heights, so calm wind stays zero and
the canopy top remains continuous. Ground wind is set to zero as in the source.
heights uses exact requested sampling points; otherwise dz sets approximate
layer spacing. For SPAC, all z_air points are preserved, including uneven grids.
"""
function wind_profile(; hc=21,pai=5.051612734794617,zref=46,uref=nothing,
    wind_speed=nothing,wind_u=nothing,wind_v=nothing,wind_height=nothing,
    dz=0.5,heights=nothing,L=Inf,fluxes=nothing,gradients=nothing,
    minWind=DEFAULT_MIN_WIND,stability_limits=:source,solver_rtol=1e-8,maxiter=100,
    outputDir=joinpath(@__DIR__,"wind_profile_results"),
    saveResults::Bool=true,makePlot::Bool=true,showPlot::Bool=true,verbose::Bool=true)

    dz = checked_number("dz",dz; positive=true)
    minWind = checked_number("minWind",minWind; nonnegative=true)
    solver_rtol = checked_number("solver_rtol",solver_rtol; positive=true)
    maxiter isa Integer && maxiter > 0 || throw(ArgumentError("maxiter must be a positive integer."))
    stability_limits in (:source,:none) || throw(ArgumentError("stability_limits must be :source or :none."))
    outputDir isa AbstractString && !isempty(outputDir) || throw(ArgumentError("outputDir must be a nonempty string."))
    speed,zwind,wind_source = measured_wind(uref,wind_speed,wind_u,wind_v,wind_height,zref)
    zwind > hc || throw(ArgumentError("wind_height must exceed hc; within-canopy measurements need a different inverse model."))
    zbound = max(zref,zwind)
    if L === :auto
        wind_source == :demo_default && throw(ArgumentError("L=:auto requires a measured wind input, not the demo wind."))
        state,diagnostics = automatic_state(hc,pai,zref,speed,zwind,zbound,
            fluxes,gradients,stability_limits,solver_rtol,maxiter)
    else
        fluxes === nothing && gradients === nothing ||
            throw(ArgumentError("Use L=:auto with thermal inputs, or omit thermal inputs for prescribed L."))
        state,L_requested,limited = prescribed_state(hc,pai,zbound,L,stability_limits)
        diagnostics = (; mode=:prescribed,status=limited ? :limited : :prescribed,
            closure_converged=missing,iterations=0,L_physical=L_requested,
            limited,inverse_L_residual=missing,thermal=nothing)
    end
    L_used = state.invL == 0 ? Inf : 1/state.invL
    ustar = speed*K_VON_KARMAN()/profile_factor(zwind,state)
    uref_used = ustar/K_VON_KARMAN()*profile_factor(zref,state)
    uCanopy = ustar/state.β
    floor_used = min(minWind,uCanopy)
    diagnostics = merge(diagnostics,(; wind_source,measured_wind=speed,
                        wind_height=zwind,minWind_effective=floor_used,
                        zeta_bound=(zbound-state.d)*state.invL))
    cfg = (; hc,pai,zref,uref=uref_used,dz,L=L_used,minWind,
             stability_limits,outputDir=String(outputDir))
    z = height_grid(hc,zref,dz,heights)
    raw = map(z) do zz
        zz == 0 ? 0.0 : zz <= hc ? uCanopy*exp((zz-hc)/state.decayLength) :
        ustar/K_VON_KARMAN()*profile_factor(zz,state)
    end
    u = [0 < zz < hc ? max(uu,floor_used) : uu for (zz,uu) in zip(z,raw)]
    all(isfinite,u) && all(x -> x >= 0,u) || error("Nonfinite or negative wind encountered.")
    isapprox(u[end],uref_used; rtol=1e-10,atol=1e-12) || error("Reference wind check failed.")
    diagnostics.limited && @warn "Stability bound is active; inspect L_physical, config.L, closure_converged, and inverse_L_residual." L_physical=diagnostics.L_physical L_used
    profile = (; Height_m=z,Wind_m_s=u,Wind_unfloored_m_s=raw,Relative_height=z./hc)
    result = (; profile,config=cfg,ustar,uCanopy,displacementHeight=state.d,
                beta=state.β,Lc=state.Lc,decayLength=state.decayLength,
                eta=state.eta,diagnostics)
    saveResults && save_results(result)
    if makePlot
        if !isdefined(@__MODULE__,:Plots)
            @eval import Plots
        end
        Base.invokelatest(draw_profile,z,u,cfg,saveResults,showPlot)
    end
    if verbose
        # @printf("Reference wind: %.6f m/s at %.3f m; u*: %.6f m/s; U(hc): %.6f m/s\n",uref_used,zref,ustar,uCanopy)
        # @printf("L used: %g m; status: %s; L physical: %g m\n",L_used,string(diagnostics.status),diagnostics.L_physical)
    end
    return result
end

function save_results(result)
    cfg,p = result.config,result.profile
    mkpath(cfg.outputDir)
    open(joinpath(cfg.outputDir,"wind_profile.csv"),"w") do io
        println(io,"Height_m,Wind_m_s,Wind_unfloored_m_s,Relative_height")
        for i in eachindex(p.Height_m)
            # @printf(io,"%.17g,%.17g,%.17g,%.17g\n",p.Height_m[i],p.Wind_m_s[i],
            #         p.Wind_unfloored_m_s[i],p.Relative_height[i])
        end
    end
    open(joinpath(cfg.outputDir,"wind_profile_diagnostics.txt"),"w") do io
        println(io,"Model: supplied Fortran RSL/v1; outer source bounds; exact neutral limit")
        for (name,value) in pairs(cfg)
            println(io,name," = ",value)
        end
        for name in (:ustar,:uCanopy,:displacementHeight,:beta,:Lc,:decayLength,:eta)
            println(io,name," = ",getproperty(result,name))
        end
        for (name,value) in pairs(result.diagnostics)
            println(io,name," = ",value)
        end
    end
    return nothing
end

function _wind_profile_from_spac(spac; uref=nothing,wind_speed=nothing,
    wind_u=nothing,wind_v=nothing,wind_height=nothing,heights=nothing,kwargs...)
    zs,z_air = spac.plant.zs,spac.info.z_air
    !isempty(zs) || throw(ArgumentError("spac.plant.zs must not be empty."))
    length(z_air) >= 2 || throw(ArgumentError("spac.info.z_air needs at least two heights."))
    all(isfinite,z_air) && all(diff(z_air) .> 0) ||
        throw(ArgumentError("spac.info.z_air must be finite and strictly increasing."))
    measured_wind(uref,wind_speed,wind_u,wind_v,wind_height,z_air[end]; allow_demo_default=false)
    return wind_profile(; hc=zs[end],pai=spac.canopy.lai+spac.canopy.sai,
        zref=z_air[end],dz=z_air[2]-z_air[1],heights=heights === nothing ? z_air : heights,
        uref,wind_speed,wind_u,wind_v,wind_height,kwargs...)
end

wind_profile(spac; kwargs...) = _wind_profile_from_spac(spac; kwargs...)

function wind_profile(cfg::NamedTuple; kwargs...)
    if hasproperty(cfg,:plant) && hasproperty(cfg,:canopy) && hasproperty(cfg,:info)
        return _wind_profile_from_spac(cfg; kwargs...)
    end
    return wind_profile(; cfg...,kwargs...)
end
