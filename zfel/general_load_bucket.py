import numpy as np

def general_load_bucket(npart,gbar,delg,iopt,Ns,coopLength,resWavelength,particle_position,s_steps,dels):
    '''
    random initialization of the beam load_bucket
    inputs:
    npart               # n-macro-particles per bucket
    gbar            # scaled detune parameter
    delg            # Gaussian energy spread in units of rho
    iopt            # 5=SASE, 4=seeded
    Ns              # N electrons per s-slice
    coopLength      # cooperation length
    resWavelength               # resonant wavelength
    particle_position           # particle information with positions in meter and gamma
    dels            # integration step in s0
    outputs:
    thet_init            # all buckets macro particles position
    gam_init             # all buckets macro particles energy
    '''
    thet_init=np.zeros((s_steps,npart))
    gam_init=np.zeros((s_steps,npart))
    if particle_position is None:
        for j in range(s_steps):
            [thet0,gam0] = load_bucket(npart,gbar,delg,iopt,Ns)     # load each bucket
            thet_init[j,:]=thet0
            gam_init[j,:]=gam0
    else:
        s_all=particle_position[:,0]
        gam_all=particle_position[:,1]
        for j in range(s_steps):
            thet_init[j,:]=(s_all[(j*npart):((j+1)*npart)]-j*dels*coopLength)*(2*np.pi)/resWavelength
            gam_init[j,:]=gam_all[(j*npart):((j+1)*npart)]

    return thet_init,gam_init


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

    print('load random bucket!!!')

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


