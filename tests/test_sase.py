import numpy as np
import matplotlib.pyplot as plt
import random
import datetime
import time
from zfel import sase1d


def test_sase():
    sase_input = dict(
    npart   = 512,                       # n-macro-particles per bucket 
    s_steps = 200,                      # n-sample points along bunch length
    z_steps = 200,                      # n-sample points along undulator
    energy  = 4313.34e6,                # electron energy [eV]
    eSpread = 0,                       # relative rms energy spread [1]
    emitN   = 1.2e-6,                    # normalized transverse emittance [m-rad]
    currentMax = 3400,                   # peak current [Ampere]
    beta = 26,                          # mean beta [meter]
    unduPeriod = 0.03,                 # undulator period [meter]
    unduK = 3.5 ,          # undulator parameter, K [1], array could taper. 
    unduL = 70,                         # length of undulator [meter]
    radWavelength=None,                 # Will calculate based on resonance condition for unduK[0]
    random_seed=31,                     # for reproducibility
    particle_position=None, #np.genfromtxt('./Inputs/particle_position.csv', delimiter=',') # or None,
    hist_rule='square-root',             # 'square-root' or 'sturges' or 'rice-rule' or 'self-design', number \
                                       #  of intervals to generate the histogram of eta value in a bucket
    iopt='sase',
    P0 = 0                            # small seed input power [W]
)
    
    output = sase1d.sase(sase_input)