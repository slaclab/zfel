import numpy as np

def load_bucket(n,gbar,delg,iopt,Ns):
    '''
    random initialization of the beam load_bucket
    inputs:
    n               # n-macro-particles per bucket
    gbar            # scaled detune parameter
    delg            # Gaussian energy spread in units of rho
    iopt            # 5=SASE, 4=seeded
    Ns              # N electrons per s-slice ??
    outputs:
    thet            # bucket macro particles position
    gam             # bucket macro particles energy
    '''
    nmax = 10000;
    if n>nmax:
        raise ValueError('increase nmax, subr load')

    gam=np.zeros(n)
    thet=np.zeros(n)
    if iopt==4:
        M=128                                               # number of particles in each beamlet
        nb= int(np.round(n/M))                              # number of beamlet via Fawley between 64 to 256 (x16=1024 to 4096)
        if M*nb!=n:
            raise ValueError('n must be a multiple of 4')
        for i in range(nb):
            #gamma=delg*np.random.randn(1)+gbar
            gamma=delg*(np.random.rand(1)-0.5)+gbar
            for j in range(M):
                gam[i*M+j]=gamma
                thet[i*M+j]=2*np.pi*(j+1)/M
    elif iopt==5:
        M=32  # number of particles in each beamlet
        nb= int(np.round(n/M) )    #number of beamlet via Fawley between 64 to 256 (x16=1024 to 4096)
        if M*nb!=n:
            raise ValueError('n must be a multiple of 4')
        effnoise = np.sqrt(3*M/(Ns/nb))    # Penman algorithm for Ns/nb >> M
        for i in range(nb):
            #gamma=delg*np.random.randn(1)+gbar
            gamma=delg*(np.random.rand(1)-0.5)+gbar
            for j in range(M):
                gam[i*M+j]=gamma
                thet[i*M+j]=2*np.pi*(j+1)/M+2*np.random.rand(1)*effnoise
    return thet,gam