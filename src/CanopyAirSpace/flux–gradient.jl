"""

    ϕ_m(ζ::FT) where {FT<:AbstractFloat}
    ϕ_m(ζ::FT, ξ::FT, c₁::FT, c₂::FT) where {FT<:AbstractFloat}

Flux–gradient relationship of Monin–Obukhov similarity theory for momentum, given
- `ζ` the stability parameter (z-d)/Lmo.
- `ξ` additional dimensionless parameter (z−d) * β / lₘ
- `c₁` coefficient for the exponential correction term (eq A16)
- `c₂` coefficient for the exponential correction term (0.5)

Equation from Bonan et al. (2018) CLM-ml (eq A10, A14, A15, A16)

"""
function ϕ_m end;

ϕ_m(ζ::FT) where {FT<:AbstractFloat} = (
    return if (ζ < 0)
        (1 - 16*ζ)^(-FT(0.25))
    else
        1 + 5*ζ
    end
);

ϕ_m(ζ::FT, ξ::FT, c₁::FT, c₂::FT) where {FT<:AbstractFloat} = (
    ϕₘ = ϕ_m(ζ);
    ϕ̂ₘ = 1 - c₁ * exp(-c₂ * ξ);

    return ϕₘ * ϕ̂ₘ
);


"""

    ϕ_c(ζ::FT) where {FT<:AbstractFloat}
    Φ_c(ζ::FT, ξ::FT, c₁::FT, c₂::FT) where {FT<:AbstractFloat}

Flux–gradient relationship of Monin–Obukhov similarity theory for heat and vapor, given
- `ζ` is the stability parameter (z-d)/Lmo
- `ξ` additional dimensionless parameter (z−d) * β / lₘ
- `c₁` coefficient for the exponential correction term (eq A19)
- `c₂` coefficient for the exponential correction term (0.5)

Equation from Bonan et al. (2018) CLM-ml (eq A11, A15, A18, A19)

"""
function ϕ_c end;

ϕ_c(ζ::FT) where {FT<:AbstractFloat} = (
    return if (ζ < 0)
        (1 - 16*ζ)^(-FT(0.5))
    else
        1 + 5*ζ
    end
);

ϕ_c(ζ::FT, ξ::FT, c₁::FT, c₂::FT) where {FT<:AbstractFloat} = (
    ϕc = ϕ_c(ζ);
    ϕ̂c = 1 - c₁ * exp(-c₂ * ξ);

    return ϕc * ϕ̂c
);



"""
eq A17
z′ is the dimensionless height within the canopy z-d
"""
function dψ̂_dt(ϕ::Function, z′::FT, Lmo::FT, lₘ::FT, β::FT, c₁::FT, c₂::FT) where {FT<:AbstractFloat}
    return ϕ(z′ / Lmo) * c₁ * exp(-c₂ * z′ * β / lₘ) / z′
end;
