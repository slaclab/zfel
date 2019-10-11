import numpy as np
import matplotlib.pyplot as plt

def general_load_bucket(npart,gbar,delg,iopt,Ns,coopLength,resWavelength,particle_position,s_steps,dels,hist_steps):
    '''
    random initialization of the beam load_bucket
    inputs:
    npart               # n-macro-particles per bucket
    gbar            # scaled detune parameter
    delg            # Gaussian energy spread in units of rho
    iopt            # 5=SASE, 4=seeded
    Ns              # N electrons per s-slice at maximum current
    coopLength      # cooperation length
    resWavelength               # resonant wavelength
    particle_position           # particle information with positions in meter and gamma
    s_steps         # n-sample points along bunch length
    dels            # integration step in s0
    hist_steps      # number of intervals to generate the histogram of gamma value in a bucket

    outputs:
    thet_init            # all buckets macro particles position
    gam_init             # all buckets macro particles energy
    '''
    if particle_position is None:
        thet_init=np.zeros((s_steps,npart))
        gam_init=np.zeros((s_steps,npart))
        for j in range(s_steps):
            [thet0,gam0] = load_bucket(npart,gbar,delg,iopt,Ns)     # load each bucket
            thet_init[j,:]=thet0
            gam_init[j,:]=gam0
    else:
        #load particle information and classify them to different intervals
        s_all=particle_position[:,0]
        #gam_all=particle_position[:,1]
        gam_all=np.random.randn(s_steps*npart)*0.01+100
        N_input=np.zeros(s_steps)
        gam_step=[[] for x in range(s_steps)]
        for k in range(s_all.shape[0]):
            location=int(s_all[k]/(dels*coopLength))
            N_input[location]+=1
            gam_step[location].append(gam_all[k])
        N_real=N_input/np.max(N_input)*Ns
        #generate theta and gamma
        thet_init=make_theta(npart,N_real)
        gam_init=make_gamma(gam_step,npart,hist_steps)

    return thet_init,gam_init




def load_bucket(n,gbar,delg,iopt,Ns):
    '''
    random initialization of the beam load_bucket
    inputs:
    n               # n-macro-particles per bucket
    gbar            # scaled detune parameter
    delg            # Gaussian energy spread in units of rho
    iopt            # 5=SASE, 4=seeded
    Ns              # N electrons per s-slice
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

def make_theta(n,N_real):
    '''
    random initialization of the beam load_bucket
    inputs:
    n               # n-macro-particles per bucket
    N_real          # real number of particles in each bucket
    outputs:
    thet_init       # all buckets macro particles position
    '''
    
    steps=N_real.shape[0]
    thet_init=np.zeros((steps,n))
    M=32  # number of particles in each beamlet
    nb= int(np.round(n/M) )    #number of beamlet via Fawley between 64 to 256 (x16=1024 to 4096)
    if M*nb!=n:
        raise ValueError('n must be a multiple of 4')
        
    for k in range(steps):
        effnoise = np.sqrt(3*M/(N_real[k]/nb))    # Penman algorithm for Ns/nb >> M
        for i in range(nb):       
            for j in range(M):
                thet_init[k,i*M+j]=2*np.pi*(j+1)/M+2*np.random.rand(1)*effnoise
    return thet_init



def make_gamma(gam_step,npart,hist_steps):
    '''
    gam_step        # input particles' gamma values in different intervals along the beam
    npart           # n-macro-particles per bucket
    hist_steps      # number of intervals to generate the histogram of gamma value in a bucket
    outputs:
    gam_init        # all buckets macro particles energy
    '''
    steps=len(gam_step)
    if np.min(gam_step)==np.max(gam_step):
        #deal with the case that all input gamma values are the same
        gam_init=np.ones((steps,npart))*np.max(gam_step)
        print('Warning, all input gamma values are the same!')
    else:
        #make histogram for gamma
        lowbound=np.min(gam_step)
        upbound=np.max(gam_step)+1e-10
        gam_hist=np.zeros((steps,hist_steps))
        for j in range(steps):
            gam_hist[j],bins=np.histogram(np.array(gam_step[j]),bins=np.linspace(lowbound, upbound, num=hist_steps+1))
        plt.figure()
        #_=plt.hist(np.array(gam_step[0]),bins=np.linspace(lowbound, upbound, num=hist_steps+1))
        #plt.title('Input gamma histogram')
        gam_hist/=np.sum(gam_hist,axis=1,keepdims=True)

        #make cdf
        gam_cdf=np.zeros(gam_hist.shape)
        gam_cdf[:,0]=gam_hist[:,0]
        for j in range(1,hist_steps):
            gam_cdf[:,j]=gam_hist[:,j]+gam_cdf[:,j-1]
        gam_cdf=np.concatenate((np.zeros((steps,1)),gam_cdf),axis=1)

        #make gamma
        x=np.random.rand(steps,npart)
        gam_init=np.zeros((steps,npart))
        for j in range(steps):
            gam_init[j]=np.interp(x[j], gam_cdf[j], bins)
        #plt.figure()
        #_=plt.hist(gam_init[0],bins=np.linspace(lowbound, upbound, num=hist_steps+1))
        #plt.title('Sampled gamma histogram')
    return gam_init


