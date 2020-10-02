import numpy as np
import matplotlib.pyplot as plt

def general_load_bucket(npart,Ns,coopLength,particle_position,s_steps,dels,hist_rule,gbar=None,delg=None,iopt=None,):
    '''
    random initialization of the beam load_bucket
    inputs:
    npart               # n-macro-particles per bucket
    gbar                # scaled detune parameter
    delg                # Gaussian energy spread in units of rho
    iopt                # 'sase' or 'seeded'
    Ns                  # N electrons per s-slice at maximum current
    coopLength          # cooperation length
    particle_position   # particle information with positions in meter and eta
    s_steps             # n-sample points along bunch length
    dels                # integration step in s0
    hist_rule           # different rules to select number of intervals to generate the histogram of eta value in a bucket

    outputs:
    thet_init           # all buckets macro particles position
    eta_init            # all buckets macro particles relative energy
    N_real              # real number of particles along the beam
    '''
    if particle_position is None:
        thet_init=np.zeros((s_steps,npart))
        eta_init=np.zeros((s_steps,npart))
        for j in range(s_steps):
            [thet0,eta0] = load_bucket(npart,gbar,delg,iopt,Ns)     # load each bucket
            thet_init[j,:]=thet0
            eta_init[j,:]=eta0
        N_real=np.ones(s_steps)
    else:
        #load particle information and classify them to different intervals
        s_all=particle_position[:,0]
        eta_all=particle_position[:,1]
        s_steps=int(np.max(s_all)/(dels*coopLength))+1 if np.max(s_all)%(dels*coopLength)!=0 else np.max(s_all)/(dels*coopLength)
        N_input=np.zeros(s_steps)
        eta_step=[[] for x in range(s_steps)]
        for k in range(s_all.shape[0]):
            location=int(s_all[k]/(dels*coopLength))
            N_input[location]+=1
            eta_step[location].append(eta_all[k])
        N_real=N_input/np.max(N_input)*Ns
        #generate theta and eta
        thet_init=np.zeros((s_steps,npart))
        eta_init=np.zeros((s_steps,npart))
        for k in range(s_steps):
            if N_real[k]==0:
                thet_init[k,:]=np.random.rand(1)*2*np.pi
                eta_init[k,:]=np.zeros(npart)
            else:
                thet_init[k,:]=make_theta(npart,N_real[k])
                eta_init[k,:]=make_eta(eta_step[k],npart,hist_rule)

    return {'thet_init':thet_init,'eta_init':eta_init,'N_real':N_real,'s_steps':s_steps}




def load_bucket(n,gbar,delg,iopt,Ns):
    '''
    random initialization of the beam load_bucket
    inputs:
    n               # n-macro-particles per bucket
    gbar            # scaled detune parameter
    delg            # Gaussian energy spread in units of rho
    iopt            # 'sase' or 'seeded'
    Ns              # N electrons per s-slice
    outputs:
    thet            # bucket macro particles position
    eta             # bucket macro particles relative energy
    '''
    nmax = 10000;
    if n>nmax:
        raise ValueError('increase nmax, subr load')

    #print('load random bucket!!!')

    eta=np.zeros(n)
    thet=np.zeros(n)
    if iopt=='seeded':
        M=128                                               # number of particles in each beamlet
        nb= int(np.round(n/M))                              # number of beamlet via Fawley between 64 to 256 (x16=1024 to 4096)
        if M*nb!=n:
            raise ValueError('n must be a multiple of 4')
        for i in range(nb):
            etaa=delg*np.random.randn(1)+gbar
            #etaa=delg*(np.random.rand(1)-0.5)+gbar
            for j in range(M):
                eta[i*M+j]=etaa
                thet[i*M+j]=2*np.pi*(j+1)/M
    elif iopt=='sase':
        M=32  # number of particles in each beamlet
        nb= int(np.round(n/M) )    #number of beamlet via Fawley between 64 to 256 (x16=1024 to 4096)
        if M*nb!=n:
            raise ValueError('n must be a multiple of 4')
        effnoise = np.sqrt(3*M/(Ns/nb))    # Penman algorithm for Ns/nb >> M
        for i in range(nb):
            etaa=delg*np.random.randn(1)+gbar
            #etaa=delg*(np.random.rand(1)-0.5)+gbar
            for j in range(M):
                eta[i*M+j]=etaa
                thet[i*M+j]=2*np.pi*(j+1)/M+2*np.random.rand(1)*effnoise
    return thet,eta

def make_theta(n,N_real_bucket):
    '''
    random initialization of a bucket's particle positions
    inputs:
    n               # n-macro-particles per bucket
    N_real_bucket   # real number of particles in a bucket
    outputs:
    thet            # macro particles position in a bucket
    '''
    
    thet=np.zeros(n)
    M=32  # number of particles in each beamlet
    nb= int(np.round(n/M) )    #number of beamlet via Fawley between 64 to 256 (x16=1024 to 4096)
    if M*nb!=n:
        raise ValueError('n must be a multiple of 4')
        
    effnoise = np.sqrt(3*M/(N_real_bucket/nb))    # Penman algorithm for Ns/nb >> M
    for i in range(nb):       
        for j in range(M):
            thet[i*M+j]=2*np.pi*(j+1)/M+2*np.random.rand(1)*effnoise
    return thet



def make_eta(eta_step_bucket,npart,hist_rule='square-root'):
    '''
    eta_step_bucket     # input particles' eta values in a bucket
    npart               # n-macro-particles per bucket
    hist_rule          # different rules to select number of intervals to generate the histogram of eta value in a bucket
    outputs:
    eta_sampled         # sampled macro particles relative energy in a bucket
    '''

    lowbound=np.min(eta_step_bucket)
    upbound=np.max(eta_step_bucket)+1e-10
    pts=len(eta_step_bucket)
    if hist_rule == 'square-root':
        hist_num=int(np.sqrt(pts))
    elif hist_rule=='sturges':
        hist_num=int(np.log2(pts))+1
    elif hist_rule=='rice-rule':
        hist_num=int(2*pts**(1/3))
    eta_hist=np.zeros(hist_num)
    eta_hist,bins=np.histogram(np.array(eta_step_bucket),bins=np.linspace(lowbound, upbound, num=hist_num+1))
        #plt.figure()
        #_=plt.hist(np.array(eta_step_bucket),bins=np.linspace(lowbound, upbound, num=hist_num+1))
        #plt.title('Input eta histogram')
    eta_hist=eta_hist/np.sum(eta_hist)

    #make cdf
    eta_cdf=np.zeros(eta_hist.shape)
    eta_cdf[0]=eta_hist[0]
    for j in range(1,hist_num):
        eta_cdf[j]=eta_hist[j]+eta_cdf[j-1]
    eta_cdf=np.concatenate((np.zeros(1),eta_cdf))
        
    #make eta
    x=np.random.rand(npart)
    eta_sampled=np.interp(x, eta_cdf, bins)
        #plt.figure()
        #_=plt.hist(eta_sampled,bins=np.linspace(lowbound, upbound, num=hist_num+1))
        #plt.title('Sampled eta histogram')
    return eta_sampled


