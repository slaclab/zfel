import numpy as np
from scipy.special import jv


def mingxie(
    sigma_x=None,  # RMS beam size
    und_lambda=None,  # Undulator period (m)
    und_k=None,  # Undulator K
    current=None,  # Beam current (A)
    gamma=None,  # Relativistic gamma
    norm_emit=None,  # Normalized emittance (m-rad)
    sigma_E=None,  # RMS energy spread (eV)
):
    """
    Calculates gain length using the Ming Xie formula.

    Design Optimization for an X-ray Free Elecgron Laser Driven by SLAC Linac
    Ming Xie, Lawrence Berkeley Laboratory, Berkeley, CA 94720, USA
    http://accelconf.web.cern.ch/AccelConf/p95/ARTICLES/TPG/TPG10.PDF

    Inputs keyword arguments:
        sigma_x      # RMS beam size
        und_lambda   # Undulator period (m)
        und_k        # Undulator K
        current      # Beam current (A)
        gamma        # Relativistic gamma
        norm_emit    # Normalized emittance (m-rad)
        sigma_E      # RMS energy spread (eV)


    Output as dict:
        gain_length        # Gain length (m)
        saturation_length  # Saturation length (m)
        saturation_power   # Saturation power (W)
        fel_wavelength     # FEL wavelength (m)
        pierce_parameter   # Pierce parameter (rho)
    """

    mec2 = 0.51099895000e6  # Electron rest mass in eV
    c_light = 299792458.0

    I = current
    IA = 17045.0
    sigma_y = sigma_x
    sigmae = sigma_E
    E = gamma * mec2

    felwave = und_lambda * (1 + und_k**2 / 2) / 2 / gamma**2

    emittance = norm_emit / gamma

    beta = sigma_x**2 / emittance
    # beta function

    # calc pierce parameter
    unduJJ = jv(0, und_k**2 / (4 + 2 * und_k**2)) - jv(
        1, und_k**2 / (4 + 2 * und_k**2)
    )
    rho = (
        (I / IA)
        * ((und_lambda * und_k * unduJJ) / (2 * np.pi)) ** 2
        / (2 * sigma_x * sigma_y)
    ) ** (1.0 / 3) / (2 * gamma)

    L1d = und_lambda / (4 * np.pi * rho * 3**0.5)
    # 1D gain length

    Lr = 4 * np.pi * sigma_x**2 / felwave
    # Rayleigh length

    yd = L1d / Lr
    # Xie's three parameters for 3D
    yr = (4 * np.pi * L1d * sigma_E) / (und_lambda * E)
    ye = (((L1d / beta) * 4 * np.pi) * emittance) / felwave
    y = (
        0.45 * yd**0.57
        + 0.55 * ye**1.6
        + 3 * yr**2
        + 0.35 * (ye**2.9) * yr**2.4
        + 51 * (yd**0.95) * yr**3
        + 5.4 * (yd**0.7) * ye**1.9
        + 1140 * ((yd**2.2) * ye**2.9) * yr**3.2
    )
    Lg = (1 + y) * L1d

    Pbeam = I * E
    Psat = (1.6 * rho * (L1d / Lg) ** 2) * Pbeam
    alpha = 1 / 9
    Pn = rho**2 * E * 1.6e-19 * c_light / felwave
    # Eq. 3, and converting eV -> J
    Lsat = Lg * np.log(Psat / Pn / alpha)  # Eq. 2

    d = {
        "gain_length": Lg,
        "saturation_length": Lsat,
        "saturation_power": Psat,
        "fel_wavelength": felwave,
        "pierce_parameter": rho,
    }
    return d
    # return Lg, Psat, felwave, rho, Lsat;
