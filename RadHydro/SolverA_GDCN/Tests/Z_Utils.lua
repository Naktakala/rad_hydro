a_const = 0.013722354852
c_const = 299.792

function SoundSpeedGiven_e_gamma(e, gamma)
    return math.sqrt(e * gamma * (gamma-1))
end

function InternalEGiven_p_gamma_rho(p, gamma, rho)
    e = p/(gamma-1)/rho
    return e
end

function InternalEGiven_T_Cv(T,Cv)
    e = Cv*T
    return e
end

function PressureGivenGammaRhoInternalE(gamma, rho, e)
    p = (gamma-1)*rho*e
    return p
end