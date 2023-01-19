import numpy as np
import scipy
from scipy import special

from zfel.particles import general_load_bucket
from zfel.fel import FEL_process_complex, final_calc

# Some constant values
alfvenCurrent = 17045.0  # Alfven current ~ 17 kA
mc2 = 0.51099906e6  # 510.99906E-3      # Electron rest mass in eV
c = 2.99792458e8  # light speed in meter
e = 1.60217733e-19  # electron charge in Coulomb
epsilon_0 = 8.85418782e-12  # electric constant
hbar = 6.582e-16  # in eV


def sase(inp_struct):
    """
    SASE 1D FEL run function
    TODO: needs updating
    Input:
    npart                       # n-macro-particles per bucket
    s_steps                     # n-sample points along bunch length
    z_steps                     # n-sample points along undulator
    energy                      # electron energy [eV]
    eSpread                     # relative rms energy spread [1]
    emitN                       # normalized transverse emittance [m-rad]
    currentMax                  # peak current [Ampere]
    beta                        # mean beta [meter]
    unduPeriod                  # undulator period [meter]
    unduK                       # undulator parameter, K [1]
    unduL                       # length of undulator [meter]
    radWavelength               # seed wavelength? [meter], used only in single-freuqency runs
    dEdz                        # rate of relative energy gain or taper [keV/m], optimal~130
    iopt                        # 'sase' or 'seeded'
    P0                          # small seed input power [W]
    random_seed                 # A random number seed. Default: None
    particle_position           # particle information with positions in meter and eta. Default: None
    hist_rule                   # different rules to select number of intervals to generate the histogram of eta value in a bucket

    Output:
    (TODO: needs updating)
    z                           # longitudinal steps along undulator
    power_z                     # power profile along undulator
    s                           # longitudinal steps along beam
    power_s                     # power profile along beam
    rho                         # FEL Pierce parameter
    detune                      # deviation from the central energy
    field                       # final output field along beam
    field_s                     # output field along beam for different z position
    gainLength                  # 1D FEL gain Length
    resWavelength               # resonant wavelength
    thet_out                    # output phase
    eta_out                     # output energy in unit of mc2
    bunching                    # bunching factor
    spectrum                    # spectrum power
    freq                        # frequency in eV
    Ns                          # real number of examples
    """

    # calculating intermediate parameters
    params = params_calc(**inp_struct)

    # Load Buckets
    bucket_params = {
        "npart": inp_struct["npart"],
        "Ns": params["Ns"],
        "coopLength": params["coopLength"],
        "particle_position": inp_struct["particle_position"],
        "s_steps": inp_struct["s_steps"],
        "dels": params["dels"],
        "hist_rule": inp_struct["hist_rule"],
        "gbar": params["gbar"],
        "delg": params["delg"],
        "iopt": inp_struct["iopt"],
    }
    bucket_data = general_load_bucket(**bucket_params)

    # Convenience
    i = inp_struct
    p = params

    # FEL process
    FEL_data = FEL_process_complex(
        i["npart"],
        i["z_steps"],
        p["kappa_1"],
        p["density"],
        p["Kai"],
        p["ku"],
        p["delt"],
        p["dels"],
        p["deta"],
        bucket_data["thet_init"],
        bucket_data["eta_init"],
        bucket_data["N_real"],
        i["s_steps"],
        E02=p["E02"],
    )

    # Finalize
    final_data = final_calc(
        FEL_data["Er"],
        FEL_data["Ei"],
        # FEL_data['eta'],
        i["s_steps"],
        i["z_steps"],
        p["kappa_1"],
        p["density"],
        p["Kai"],
        p["Pbeam"],
        p["delt"],
        p["dels"],
    )

    # Collect output
    output = FEL_data
    output.update(final_data)
    output["params"] = params

    # Extra (put somewhere else)
    s = (
        np.arange(1, i["s_steps"] + 1) * p["dels"] * p["coopLength"]
    )  # longitundinal steps along beam in m
    z = (
        np.arange(1, i["z_steps"] + 1) * p["delt"]
    )  # longitundinal steps along undulator in meter
    bunchLength = s[-1]  # beam length in meter
    bunch_steps = np.round(
        bunchLength / p["delt"] / p["coopLength"]
    )  # rms (Gaussian) or half width (flattop) bunch length in s_step
    output["s"] = s
    output["z"] = z
    output["bunchLength"] = bunchLength
    output["bunch_steps"] = bunch_steps

    omega = hbar * 2.0 * np.pi / (p["resWavelength"] / c)
    df = hbar * 2.0 * np.pi * 1 / (bunchLength / c)
    output["freq"] = np.linspace(
        omega - i["s_steps"] / 2 * df, omega + i["s_steps"] / 2 * df, i["s_steps"]
    )

    return output


def params_calc(
    *,  # Require kwargs explicitly to avoid typos
    npart=512,
    s_steps=200,
    z_steps=200,
    energy=4313.34e6,
    eSpread=0,
    emitN=1.2e-6,
    currentMax=3400,
    beta=26,
    unduPeriod=0.03,
    unduK=3.5,
    unduL=70,
    iopt="sase",
    P0=0,
    random_seed=None,
    particle_position=None,
    radWavelength=None,
    hist_rule="square-root"
):
    """
    calculating intermediate parameters
    """
    # whether to use constant random seed for reproducibility
    if random_seed is not None:
        np.random.seed(random_seed)

    # Check if unduK is array. Otherwise, fill it out.
    if not isinstance(unduK, np.ndarray):
        unduK = np.full(z_steps, unduK)

    unduJJ = scipy.special.jv(0, unduK**2 / (4 + 2 * unduK**2)) - scipy.special.jv(
        1, unduK**2 / (4 + 2 * unduK**2)
    )  # undulator JJ
    gamma0 = energy / mc2  # central energy of the beam in unit of mc2
    sigmaX2 = (
        emitN * beta / gamma0
    )  # rms transverse size, divergence of the electron beam

    # Needed for FEL_process
    kappa_1 = (
        e * unduK * unduJJ / 4 / epsilon_0 / gamma0
    )  # Eq. 4.10 in Kim, Huang, Lindberg (2017)
    Kai = e * unduK * unduJJ / (2 * gamma0**2 * mc2 * e)  # Ibid.
    density = currentMax / (e * c * 2 * np.pi * sigmaX2)
    ku = 2 * np.pi / unduPeriod

    rho = (0.5 / gamma0) * (
        (currentMax / alfvenCurrent)
        * (unduPeriod * unduK * unduJJ / (2 * np.pi)) ** 2
        / (2 * sigmaX2)
    ) ** (
        1 / 3
    )  # FEL Pierce parameter

    resWavelength = (
        unduPeriod * (1 + unduK[0] ** 2 / 2.0) / (2 * gamma0**2)
    )  # resonant wavelength

    if radWavelength is None:
        radWavelength = resWavelength

    Pbeam = energy * currentMax  # beam power [W]
    coopLength = resWavelength / unduPeriod  # cooperation length
    # cs0  = bunchLength/coopLength                          # bunch length in units of cooperation length
    z0 = unduL  # wiggler length
    delt = z0 / z_steps  # integration step in z0 ~ 0.1 gain length
    dels = delt  # integration step in s0 must be same as in z0
    E02 = density * kappa_1[0] * P0 / Pbeam / Kai[0]  # scaled input power
    gbar = resWavelength / radWavelength - 1.0  # scaled detune parameter

    delg = eSpread  # Gaussian energy spread in units of rho
    Ns = (
        currentMax * unduL / unduPeriod / z_steps * resWavelength / c / e
    )  # N electrons per s-slice [ ]

    deta = np.sqrt((1 + 0.5 * unduK[0] ** 2) / (1 + 0.5 * unduK**2)) - 1

    params = {
        "unduJJ": unduJJ,
        "gamma0": gamma0,
        "sigmaX2": sigmaX2,
        "kappa_1": kappa_1,
        "density": density,
        "Kai": Kai,
        "ku": ku,
        "resWavelength": resWavelength,
        "Pbeam": Pbeam,
        "coopLength": coopLength,
        "z0": z0,
        "delt": delt,
        "dels": dels,
        "E02": E02,
        "gbar": gbar,
        "delg": delg,
        "Ns": Ns,
        "deta": deta,
        "rho": rho,
    }

    return params
