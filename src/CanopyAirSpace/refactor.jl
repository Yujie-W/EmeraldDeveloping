
function canopy_state(hc::FT, pai::FT, invL::FT, Lmo::FT) where {FT<:AbstractFloat}
    # drag coefficient is set to be 0.25
    # Bonan et al. (2018) CLM-ml
    #     - Lc  : eq A27
    #     - βₙ  : eq A31 (max at 0.35)
    #     - β   : eq A24 and A23 (stable and unstable)
    #     - hd  : eq A33 (hc - d)
    #     - Sc1 : eq A25
    #     - Sc  : eq A34
    #     - Φm  : eq A10
    #     - Φc  : eq A11
    #     - c₁m : eq A16
    #     - c₁c : eq A19
    #     - lₘ  : eq A26
    #     - η   : page 1487  (max set to 20)
    vk  = K_VON_KARMAN(FT);
    Lc  = hc / pai / FT(0.25);
    βₙ  = β_neutral(hc, pai);
    β   = (Lmo > 0) ? β_stable(Lc, Lmo, βₙ) : β_unstable(Lc, Lmo, βₙ);
    hd  = min(hc, β^2 * Lc * (1 - exp(-FT(0.25) * pai / β^2)));
    d   = hc - hd;
    Sc0 = FT(1);
    Sc1 = FT(0.5) + FT(0.3) * tanh(2*Lc/Lmo);
    Sc  = (1 - βₙ/FT(0.35)) * Sc0 + βₙ/FT(0.35) * Sc1;
    ϕm  = ϕ_m(hd / Lmo);
    ϕc  = ϕ_c(hd / Lmo);
    c₂  = FT(0.5);
    c₁m = (1 - vk / (2 * β * ϕm)) * exp(c₂/2);
    c₁c = (1 - vk / (2 * β * ϕc)) * exp(c₂/2);
    lₘ  = 2 * β^3 * Lc;
    η   = min(20, hc * β / lₘ);

    return (; hc, pai, invL, Lc, βₙ, β, hd, d, Sc, c₁m, c₁c, η, decayLength=hc/η)
end
