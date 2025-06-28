using QuadGK
using NLsolve
using Printf

# Constants
const integration_accuracy = 1e-3
const integration_accuracy_a = 1e-6
const numb = 20000

const Nc = 3.0
const c2 = 0.5517
const c1 = 0.830185185
const fpi = 92.0

const m_l = 3.5
const m_h = 95.0

# Global variables (mutable)
qr_p = 0.011091
qi_p = 0.0
sl_p = 1.6572
sh_p = 40.361468

T = 0.0
mu = 0.0
mua = 0.0
T0 = 270.0

g = 4.0
Y = g  # Synonym
M = Y * fpi / 2.0

# Parameters
h_u = 121.716^3
h_s = 384.424^3
c = 4557.82
lambda = 18.2552 + 0.0396 * Y^4
m2 = 537.606^2 - 11.2915^2 * Y^4

# Helper functions
x2(x) = x^2
x3(x) = x^3

function Usigma(sl, sh, qr, qi)
    out = m2 * (2.0 * x2(sl) + x2(sh)) +
          lambda * (2.0 * sl^4 + sh^4) -
          2.0 * c * sl * sl * sh -
          2.0 * sl * h_u -
          1.0 * sh * h_s
    return out
end

function dUsigma_dsl(sl, sh, qr, qi)
    out = m2 * 4 * sl +
          lambda * 8 * sl^3 -
          c * 4 * sl * sh -
          2 * h_u
    return out
end

function dUsigma_dsldsl(sl, sh)
    out = m2 * 4 +
          lambda * 24 * sl^2 -
          c * 4 * sh
    return out
end

function dUsigma_dsh(sl, sh, qr, qi)
    out = m2 * 2 * sh +
          lambda * 4 * sh^3 -
          c * 2 * x2(sl) -
          1 * h_s
    return out
end

function dUsigma_dshdsh(sl, sh)
    out = m2 * 2 +
          lambda * 12 * sh^2
    return out
end

function dUsigma_dshdsl(sl, sh)
    return -c * 4 * sl
end

function dUsigma_dsldqi(qr, qi)
    return 0.0
end

function dUsigma_dshdqi(qr, qi)
    return 0.0
end

function dUsigma_dsldqr(qr, qi)
    return 0.0
end

function dUsigma_dshdqr(qr, qi)
    return 0.0
end

function d1func()
    return 2.0 * π^2 / 15.0 * c1 * (T0/T)^2
end

function d2func()
    return 2.0 * π^2 / 3.0 * (1 - c2 * (T0/T)^2)
end

function gluonpotential_dqr(qr, qi, sl, sh)
    d1 = d1func()
    d2 = d2func()
    glue = 8 * (d1 * (-1 + 3*qr) - 3*d2 * (-1 + 2*qr) * (9*qi^2 + qr - 3*qr^2))
    return glue
end

function gluonpotential_dqi(qr, qi, sl, sh)
    d1 = d1func()
    d2 = d2func()
    glue = -72 * qi * (d1 + d2 * (1 - 18*qi^2 - 6*qr + 6*qr^2))
    return glue
end

function gluonpotential_dqrdqr(qr, qi, sl, sh)
    d1 = d1func()
    d2 = d2func()
    glue = 24 * (d1 + d2 * (1 - 18*qi^2 - 10*qr + 18*qr^2))
    return glue
end

function gluonpotential_dqidqi(qr, qi, sl, sh)
    d1 = d1func()
    d2 = d2func()
    glue = -72 * (d1 + d2 * (1 - 54*qi^2 - 6*qr + 6*qr^2))
    return glue
end

# Thermal integrals
function pressure(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        p^2 * (log(((1 + exp(4*π*qi)*f) * (exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr))) / exp(4*π*qi)) +
               log(((exp(4*π*qi) + fa) * (1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr))) / exp(4*π*qi)))
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result
end

function nb(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        p^2 * (((f * (exp(8*π*qi) + 2*f + 3*exp(4*π*qi)*f^2 + 2*exp(2*π*qi)*(1 + 2*exp(4*π*qi)*f)*cos(2*π*qr))) / 
                ((1 + exp(4*π*qi)*f) * (exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr))) -
                (fa * (1 + 2*exp(8*π*qi)*fa + 3*exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*(exp(4*π*qi) + 2*fa)*cos(2*π*qr))) / 
                ((exp(4*π*qi) + fa) * (1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr)))) / T)
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result
end

function EoM_s(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        p^2 * (m * (-((f * (exp(8*π*qi) + 2*f + 3*exp(4*π*qi)*f^2 + 2*exp(2*π*qi)*(1 + 2*exp(4*π*qi)*f)*cos(2*π*qr))) / 
                     ((1 + exp(4*π*qi)*f) * (exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr)))) -
                    (fa * (1 + 2*exp(8*π*qi)*fa + 3*exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*(exp(4*π*qi) + 2*fa)*cos(2*π*qr))) / 
                    ((exp(4*π*qi) + fa) * (1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr)))) / (Eq * T))
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result * g
end

function EoM_qr(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        p^2 * ((-4 * exp(2*π*qi) * π * ((f + exp(4*π*qi)*fa) * (1 + f*fa) + 4*exp(2*π*qi)*f*fa*cos(2*π*qr)) * sin(2*π*qr)) /
               ((exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr)) * (1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr))))
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result
end

function EoM_qi(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        p^2 * (4 * π * ((f * (exp(8*π*qi) - f + exp(2*π*qi)*(-1 + exp(4*π*qi)*f)*cos(2*π*qr))) /
                        ((1 + exp(4*π*qi)*f) * (exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr))) +
                        (fa * (-1 + exp(8*π*qi)*fa + exp(2*π*qi)*(exp(4*π*qi) - fa)*cos(2*π*qr))) /
                        ((exp(4*π*qi) + fa) * (1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr)))))
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result
end

function dPq_dsds(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        # Simplified version of the complex expression
        p^2 * (-(((Eq*f^2*m^2*(exp(8*π*qi) + 2*f + 3*exp(4*π*qi)*f^2 + 2*exp(2*π*qi)*(1 + 2*exp(4*π*qi)*f)*cos(2*π*qr))^2) /
                ((1 + exp(4*π*qi)*f)^2 * (exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr))^2) +
                (Eq*fa^2*m^2*(1 + 2*exp(8*π*qi)*fa + 3*exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*(exp(4*π*qi) + 2*fa)*cos(2*π*qr))^2) /
                ((exp(4*π*qi) + fa)^2 * (1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr))^2))) / (Eq^3 * T^2))
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result * x2(g)
end

function dPq_dqrds(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        p^2 * ((-4*exp(2*π*qi)*(-1 + f*fa)*m*π*(8*exp(6*π*qi)*f*fa*(1 + f*fa)*cos(2*π*qr) + 
               (f + exp(4*π*qi)*fa)*(-f^2 - exp(8*π*qi)*fa^2 + exp(4*π*qi)*(1 + 6*f*fa + f^2*fa^2) + 
               2*exp(4*π*qi)*f*fa*cos(4*π*qr)))*sin(2*π*qr)) /
               (Eq*T*(exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr))^2 * 
               (1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr))^2))
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result * g
end

function dPq_dqids(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        out = p^2 * ((4*m*π*(-((f*(exp(8*π*qi) - 2*f + exp(2*π*qi)*(-1 + 2*exp(4*π*qi)*f)*cos(2*π*qr))) /
                             ((1 + exp(4*π*qi)*f)*(exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr)))) +
                            (f^2*(exp(8*π*qi) - f + exp(2*π*qi)*(-1 + exp(4*π*qi)*f)*cos(2*π*qr))*
                             (exp(8*π*qi) + 2*f + 3*exp(4*π*qi)*f^2 + 2*exp(2*π*qi)*(1 + 2*exp(4*π*qi)*f)*cos(2*π*qr))) /
                            ((1 + exp(4*π*qi)*f)^2*(exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr))^2) +
                            (fa*(1 - 2*exp(8*π*qi)*fa - exp(2*π*qi)*(exp(4*π*qi) - 2*fa)*cos(2*π*qr))) /
                            ((exp(4*π*qi) + fa)*(1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr))) +
                            (fa^2*(-1 + exp(8*π*qi)*fa + exp(2*π*qi)*(exp(4*π*qi) - fa)*cos(2*π*qr))*
                             (1 + 2*exp(8*π*qi)*fa + 3*exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*(exp(4*π*qi) + 2*fa)*cos(2*π*qr))) /
                            ((exp(4*π*qi) + fa)^2*(1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr))^2))) / (Eq*T))
        
        abs(out) < 1e-10 ? 0.0 : out
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-5)
    return result * g
end

function dPq_dsdsds(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        # Simplified version - actual expression is very complex
        p^2 * ((2*((-3*f^3*m)/(Eq*T) - (f*m*(exp(4*π*qi) + (2*cos(2*π*qr))/exp(2*π*qi)))/(Eq*T) - 
               (2*f^2*m*(exp(-4*π*qi) + 2*exp(2*π*qi)*cos(2*π*qr)))/(Eq*T))^3) /
               (1 + f^3 + f*(exp(4*π*qi) + (2*cos(2*π*qr))/exp(2*π*qi)) + 
               f^2*(exp(-4*π*qi) + 2*exp(2*π*qi)*cos(2*π*qr)))^3)
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result * x2(g) * g
end

function dPq_dsdsdqr(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        # Simplified version
        p^2 * ((2*((-3*f^3*m)/(Eq*T) - (f*m*(exp(4*π*qi) + (2*cos(2*π*qr))/exp(2*π*qi)))/(Eq*T) - 
               (2*f^2*m*(exp(-4*π*qi) + 2*exp(2*π*qi)*cos(2*π*qr)))/(Eq*T))^2 * 
               ((-4*f*π*sin(2*π*qr))/exp(2*π*qi) - 4*exp(2*π*qi)*f^2*π*sin(2*π*qr))) /
               (1 + f^3 + f*(exp(4*π*qi) + (2*cos(2*π*qr))/exp(2*π*qi)) + 
               f^2*(exp(-4*π*qi) + 2*exp(2*π*qi)*cos(2*π*qr)))^3)
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result * x2(g)
end

function dPq_dsdqrdqr(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        p^2 * (((8*f*m*π^2*cos(2*π*qr))/(exp(2*π*qi)*Eq*T) + (16*exp(2*π*qi)*f^2*m*π^2*cos(2*π*qr))/(Eq*T)) /
               (1 + f^3 + f*(exp(4*π*qi) + (2*cos(2*π*qr))/exp(2*π*qi)) + 
               f^2*(exp(-4*π*qi) + 2*exp(2*π*qi)*cos(2*π*qr))))
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result * g
end

function dPq_dqrdqr(qr, qi, s)
    func(p) = begin
        m = g * s
        Eq = sqrt(p^2 + x2(m))
        f = exp(-(Eq - mu) / T)
        fa = exp(-(Eq - mua) / T)
        
        p^2 * (8*exp(2*π*qi)*π^2*(cos(2*π*qr)*(-(f/(exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr))) - 
               fa/(1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr))) + 
               2*exp(2*π*qi)*(-(f^2/(exp(4*π*qi) + f^2 + 2*exp(2*π*qi)*f*cos(2*π*qr))^2) - 
               fa^2/(1 + exp(4*π*qi)*fa^2 + 2*exp(2*π*qi)*fa*cos(2*π*qr))^2)*sin(2*π*qr)^2))
    end
    
    result, _ = quadgk(func, 0.0, 10000.0, rtol=1e-6)
    return result
end

function dOmegadsldsl(qr, qi, sl, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 2 * T
    quark = factor_light * dPq_dsds(qr, qi, sl)
    meson = dUsigma_dsldsl(sl, sh)
    vacuum = -2.0 * Nc / (8 * x2(π)) * g^4 * x2(sl) * (7 + 12 * log(g * sl / M))
    
    SB = -m_l / Y * factor_light * dPq_dsdsds(qr, qi, sl)
    return quark + meson + vacuum + SB
end

function dOmegadshdsh(qr, qi, sl, sh)
    factor_heavy = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    meson = dUsigma_dshdsh(sl, sh)
    quark = factor_heavy * dPq_dsds(qr, qi, sh)
    vacuum = -Nc / (8 * x2(π)) * g^4 * x2(sh) * (7 + 12 * log(g * sh / M))
    SB = -m_h / Y * factor_heavy * dPq_dsdsds(qr, qi, sh)
    return quark + meson + vacuum + SB
end

function dOmegadshdsl(qr, qi, sl, sh)
    meson = dUsigma_dshdsl(sl, sh)
    return meson
end

function dOmegadsldqr(qr, qi, sl, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 2 * T
    quark = factor_light * dPq_dqrds(qr, qi, sl)
    meson = dUsigma_dsldqr(qr, qi)
    SB = -m_l / Y * factor_light * dPq_dsdsdqr(qr, qi, sl)
    return quark + meson + SB
end

function dOmegadshdqr(qr, qi, sl, sh)
    factor_heavy = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    meson = dUsigma_dshdqr(qr, qi)
    quark = factor_heavy * dPq_dqrds(qr, qi, sh)
    SB = -m_h / Y * factor_heavy * dPq_dsdsdqr(qr, qi, sh)
    return quark + meson + SB
end

function dOmegadqrdqr(qr, qi, sl, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 2 * T
    factor_heavy = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    gluon = gluonpotential_dqrdqr(qr, qi, sl, sh)
    quark = factor_heavy * dPq_dqrdqr(qr, qi, sh) + factor_light * dPq_dqrdqr(qr, qi, sl)
    SB = -m_l / Y * factor_light * dPq_dsdqrdqr(qr, qi, sl) - 
         m_h / Y * factor_heavy * dPq_dsdqrdqr(qr, qi, sh)
    return quark + gluon * x2(T) * x2(T) + SB
end

function Interface!(F, x)
    sl = abs(x[1])
    sh = abs(x[2])
    qr = x[3]
    qi = x[4]
    
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 2 * T
    factor_heavy = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    
    F[1] = factor_light * EoM_s(qr, qi, sl) / T^4 -
           (g^4 * Nc * 2 * sl^3 * (1 + 4 * log(g * sl / M))) / (8.0 * π^2) / T^4 +
           dUsigma_dsl(sl, sh, qr, qi) / T^4 -
           m_l / Y * factor_light * dPq_dsds(qr, qi, sl) / T^4
    
    F[2] = factor_light * EoM_qr(qr, qi, sl) / T^4 +
           factor_heavy * EoM_qr(qr, qi, sh) / T^4 +
           gluonpotential_dqr(qr, qi, sl, sh) -
           m_l / Y * factor_light * dPq_dqrds(qr, qi, sl) / T^4 -
           m_h / Y * factor_heavy * dPq_dqrds(qr, qi, sh) / T^4
    
    F[3] = factor_light * EoM_qi(qr, qi, sl) / T^4 +
           factor_heavy * EoM_qi(qr, qi, sh) / T^4 +
           gluonpotential_dqi(qr, qi, sl, sh) -
           m_l / Y * factor_light * dPq_dqids(qr, qi, sl) / T^4 -
           m_h / Y * factor_heavy * dPq_dqids(qr, qi, sh) / T^4
    
    F[4] = factor_heavy * EoM_s(qr, qi, sh) / T^4 -
           (g^4 * Nc * 1 * sh^3 * (1 + 4 * log(g * sh / M))) / (8.0 * π^2) / T^4 +
           dUsigma_dsh(sl, sh, qr, qi) / T^4 -
           m_h / Y * factor_heavy * dPq_dsds(qr, qi, sh) / T^4
    
    return F
end

function pressure_final(qr, qi, sl, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 2 * T
    factor_heavy = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    norm = 1.0 / T^4
    Us = Usigma(sl, sh, qr, qi) * norm
    thermal = factor_light * pressure(qr, qi, sl) * norm + factor_heavy * pressure(qr, qi, sh) * norm
    m_light = g * sl
    m_heavy = g * sh
    vacuum_light = -Nc * 2 / (8 * π^2) * m_light^4 * log(m_light / M) * norm
    vacuum_heavy = -Nc * 1 / (8 * π^2) * m_heavy^4 * log(m_heavy / M) * norm
    
    d1 = d1func()
    d2 = d2func()
    Ug = (-4*d1*(9*qi^2 + (2 - 3*qr)*qr) + 4*d2*(81*qi^4 - 9*qi^2*(1 - 6*qr + 6*qr^2) + qr^2*(3 - 10*qr + 9*qr^2)))
    
    return thermal + vacuum_light + vacuum_heavy + Us + Ug -
           m_l/Y * factor_light * norm * EoM_s(qr, qi, sl) -
           m_h/Y * factor_heavy * norm * EoM_s(qr, qi, sh)
end

function nb_final(qr, qi, sl, sh)
    factor_light = 1.0 / (2.0 * x2(π)) * 2 * 2
    factor_heavy = 1.0 / (2.0 * x2(π)) * 2 * 1
    norm = 1.0 / T^3
    thermal = factor_light * nb(qr, qi, sl) * norm + factor_heavy * nb(qr, qi, sh) * norm
    return thermal
end

function pressure_final_vacuum()
    norm = 1.0 / T^4
    sl = 46.0029
    sh = 76.0827
    qr = 0.3317
    qi = 0.0
    Us = Usigma(sl, sh, qr, qi) * norm
    m_light = g * sl
    m_heavy = g * sh
    vacuum_light = -Nc * 2 / (8 * π^2) * m_light^4 * log(m_light / M) * norm
    vacuum_heavy = -Nc * 1 / (8 * π^2) * m_heavy^4 * log(m_heavy / M) * norm
    d1 = d1func()
    d2 = d2func()
    Ug = (-4*d1*(9*qi^2 + (2 - 3*qr)*qr) + 4*d2*(81*qi^4 - 9*qi^2*(1 - 6*qr + 6*qr^2) + qr^2*(3 - 10*qr + 9*qr^2)))
    return vacuum_light + vacuum_heavy + Us + Ug
end

function SigmaMass(qr, qi, sl, sh)
    out = (dOmegadsldsl(qr, qi, sl, sh) + dOmegadshdsh(qr, qi, sl, sh) + 2.0 * dOmegadshdsl(qr, qi, sl, sh)) / 6.0
    return out / fpi / fpi
end

function PionMass(qr, qi, sl, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    H_u = h_u + m_l / Y * factor_light * dPq_dsds(qr, qi, sl)
    out = H_u / (2.0 * sl)
    return out / fpi / fpi
end

function KaonMass(qr, qi, sl, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    factor_heavy = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    
    H_u = h_u + m_l / Y * factor_light * dPq_dsds(qr, qi, sl)
    H_s = h_s + m_h / Y * factor_heavy * dPq_dsds(qr, qi, sh)
    out = (H_u + H_s) / (2.0 * (sl + sh))
    return out / fpi / fpi
end

function KappaMass(qr, qi, sl, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    factor_heavy = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    
    H_u = h_u + m_l / Y * factor_light * dPq_dsds(qr, qi, sl)
    H_s = h_s + m_h / Y * factor_heavy * dPq_dsds(qr, qi, sh)
    
    out = (H_u - H_s) / (2.0 * (sl - sh))
    return out / fpi / fpi
end

function EtaMass(qr, qi, sl, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    factor_heavy = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    
    H_u = h_u + m_l / Y * factor_light * dPq_dsds(qr, qi, sl)
    H_s = h_s + m_h / Y * factor_heavy * dPq_dsds(qr, qi, sh)
    
    sqrtdiff = sqrt((0.5*H_u/sl - 0.5*H_s/sh + 2*c*sh*(1.0 - 0.5*(sl/sh)^2))^2 + 8*c^2*sl^2)
    sum = 0.5*H_u/sl + 0.5*H_s/sh + 2*c*sh*(1.0 + 0.5*(sl/sh)^2)
    out = 0.5 * (sum - sqrtdiff)
    return out / fpi / fpi
end

function EtaPrimeMass(qr, qi, sl, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    factor_heavy = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    
    H_u = h_u + m_l / Y * factor_light * dPq_dsds(qr, qi, sl)
    H_s = h_s + m_h / Y * factor_heavy * dPq_dsds(qr, qi, sh)
    
    sqrtdiff = sqrt((0.5*H_u/sl - 0.5*H_s/sh + 2*c*sh*(1.0 - 0.5*(sl/sh)^2))^2 + 8*c^2*sl^2)
    sum = 0.5*H_u/sl + 0.5*H_s/sh + 2*c*sh*(1.0 + 0.5*(sl/sh)^2)
    out = 0.5 * (sum + sqrtdiff)
    return out / fpi / fpi
end

function DU2ds(qr, qi, sl_i, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 2 * T
    step = sl_i * 0.01
    
    sl = sl_i + step
    up = factor_light * EoM_s(qr, qi, sl) -
         (g^4 * Nc * 2 * sl^3 * (1 + 4 * log(g * sl / M))) / (8.0 * π^2) +
         dUsigma_dsl(sl, sh, qr, qi) -
         m_l / Y * factor_light * dPq_dsds(qr, qi, sl)
    
    sl = sl_i - step
    down = factor_light * EoM_s(qr, qi, sl) -
           (g^4 * Nc * 2 * sl^3 * (1 + 4 * log(g * sl / M))) / (8.0 * π^2) +
           dUsigma_dsl(sl, sh, qr, qi) -
           m_l / Y * factor_light * dPq_dsds(qr, qi, sl)
    
    return 0.5 * (up - down) / step
end

function DU2dh(qr, qi, sl, sh_i)
    factor_heavy = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    step = sh_i * 0.01
    
    sh = sh_i + step
    up = factor_heavy * EoM_s(qr, qi, sh) -
         (g^4 * Nc * 1 * sh^3 * (1 + 4 * log(g * sh / M))) / (8.0 * π^2) +
         dUsigma_dsh(sl, sh, qr, qi) -
         m_h / Y * factor_heavy * dPq_dsds(qr, qi, sh)
    
    sh = sh_i - step
    down = factor_heavy * EoM_s(qr, qi, sh) -
           (g^4 * Nc * 1 * sh^3 * (1 + 4 * log(g * sh / M))) / (8.0 * π^2) +
           dUsigma_dsh(sl, sh, qr, qi) -
           m_h / Y * factor_heavy * dPq_dsds(qr, qi, sh)
    
    return 0.5 * (up - down) / step
end

function DU2dhds(qr, qi, sl, sh)
    return -4.0 * c * sl
end

function a0Mass(qr, qi, sl_i, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 2 * T
    step = sl_i * 0.01
    
    sl = sl_i + step
    up = factor_light * EoM_s(qr, qi, sl) -
         (g^4 * Nc * 2 * sl^3 * (1 + 4 * log(g * sl / M))) / (8.0 * π^2) -
         m_l / Y * factor_light * dPq_dsds(qr, qi, sl)
    
    sl = sl_i - step
    down = factor_light * EoM_s(qr, qi, sl) -
           (g^4 * Nc * 2 * sl^3 * (1 + 4 * log(g * sl / M))) / (8.0 * π^2) -
           m_l / Y * factor_light * dPq_dsds(qr, qi, sl)
    
    out = 0.25 * 0.5 * (up - down) / step + m2 + 6.0 * lambda * sl_i^2 + c * sh
    return out / fpi / fpi
end

function chis(qr, qi, sl, sh)
    factor_light = -1.0 / (2.0 * x2(π)) * 2 * 1 * T
    
    pipi = (factor_light * EoM_s(qr, qi, sl) / (2.0 * sl) - 
            (g^4 * Nc * sl^2 * (1 + 4 * log(g * sl / M))) / (16.0 * π^2))
    a0a0 = (factor_light * dPq_dsds(qr, qi, sl) - 
            (g^4 * Nc * sl^2 * (7 + 12 * log(g * sl / M))) / (16.0 * π^2))
    
    pi_0 = -c * sh + 2 * lambda * sl^2 + m2
    a0_0 = c * sh + 6 * lambda * sl^2 + m2
    
    piressumed = 1.0 / (PionMass(qr, qi, sl, sh) * fpi * fpi)
    a0ressumed = 1.0 / (a0Mass(qr, qi, sl, sh) * fpi * fpi)
    
    chi_pipi = pipi * pi_0 * piressumed
    chi_a0a0 = a0a0 * a0_0 * a0ressumed
    
    return chi_pipi, chi_a0a0
end

function solver!(x_init)
    global qr_p, qi_p, sl_p, sh_p
    
    # Define the function for NLsolve
    f!(F, x) = Interface!(F, x)
    
    # Solve the system
    result = nlsolve(f!, x_init, method = :trust_region, ftol = 1e-12, iterations = 2000)
    
    if converged(result)
        sl_p = abs(result.zero[1])
        sh_p = abs(result.zero[2])
        qr_p = result.zero[3]
        qi_p = result.zero[4]
        return true
    else
        println("Solver did not converge")
        return false
    end
end

function main()
    global T, mu, mua, Y, g, lambda, m2, M
    global qr_p, qi_p, sl_p, sh_p
    
    # For simplicity
    Y = 5.0
    g = Y
    lambda = 18.255188335605816 + 0.03962588514828971 * Y^4
    m2 = 537.606^2 - 11.2915^2 * Y^4
    M = Y * fpi / 2.0
    
    # Open output file
    pSusc = open("chi_pipi.dat", "w")
    
    mu = 0.0
    mua = -mu
    
    for T_val in 500:-1:50
        global T = T_val
        
        if T > 3000
            break
        end
        
        mu = T * 0.0
        mua = -mu
        
        # Initial guess
        x_init = [sl_p, sh_p, qr_p, qi_p]
        
        # Solve the system
        if solver!(x_init)
            println("T = $T, sl = $sl_p, sh = $sh_p, qr = $qr_p, qi = $qi_p")
            
            chi_pipi, chi_a0a0 = chis(qr_p, qi_p, sl_p, sh_p)
            
            kappa = g^4 * 3 / 16 / π^2
            vacuum_diff = 6.0 * kappa * sl_p^2 + 8.0 * kappa * sl_p^2 * log(g * sl_p / M)
            
            @printf(pSusc, "%.25f %.25f %.25f %.25f\n", T, chi_pipi, chi_a0a0, vacuum_diff)
            flush(pSusc)
        end
    end
    
    close(pSusc)
end

# Run the main function
main()