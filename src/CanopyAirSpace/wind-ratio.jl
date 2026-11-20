# cache solver
const METHOD_NB::NewtonBisectionMethod{Float64} = NewtonBisectionMethod{Float64}();
const SOLUTION_TOL::SolutionTolerance{Float64} = SolutionTolerance{Float64}(1e-6, 50);


"""

    beta_neutral(hc::FT, pai::FT; z₀_m::FT = FT(0.01)) where {FT<:AbstractFloat}

Return the ratio of friction velocity to wind speed, given
- `hc` Canopy height
- `pai` Plant area index
- `z₀_m` Roughness length for momentum of the underlying ground surface

"""
function β_neutral(hc::FT, pai::FT; z₀_m::FT = FT(0.01)) where {FT<:AbstractFloat}
    cᵦ = (K_VON_KARMAN(FT) / log(1 + hc / z₀_m)) ^ 2;
    βₙ = sqrt(cᵦ + FT(0.3) * pai);

    return min(βₙ, FT(0.35))
end;


"""

    β_stable(Lc::FT, Lmo::FT, βₙ::FT) where {FT<:AbstractFloat}

Return the ratio of friction velocity to wind speed under stable conditions, given
- `Lc` Canopy length scale
- `Lmo` Monin-Obukhov length
- `βₙ` Neutral beta ratio

"""
function β_stable(Lc::FT, Lmo::FT, βₙ::FT) where {FT<:AbstractFloat}
    # solve 5*Lc/Lmo * β^3 + β - βₙ = 0 for β in [0, βₙ]
    METHOD_NB.x_min = 0;
    METHOD_NB.x_max = βₙ;
    METHOD_NB.x_ini = βₙ / 2;
    tarf(x) = 5*Lc/Lmo * x^3 + x - βₙ;
    β = find_zero(tarf, METHOD_NB, SOLUTION_TOL);

    return clamp(β, FT(0.2), FT(0.5))
end;


"""

    β_unstable(Lc::FT, Lmo::FT, βₙ::FT) where {FT<:AbstractFloat}

Return the ratio of friction velocity to wind speed under unstable conditions, given
- `Lc` Canopy length scale
- `Lmo` Monin-Obukhov length
- `βₙ` Neutral beta ratio

"""
function β_unstable(Lc::FT, Lmo::FT, βₙ::FT) where {FT<:AbstractFloat}
    # solve (β²)^2  + 16*Lc/Lmo*βₙ⁴ * β² - βₙ⁴ = 0 for β in [0, βₙ]
    a = FT(1);
    b = 16 * Lc / Lmo * βₙ^4;
    c = - βₙ^4;
    β² = upper_quadratic(a, b, c);

    return clamp(sqrt(β²), FT(0.2), FT(0.5))
end;
