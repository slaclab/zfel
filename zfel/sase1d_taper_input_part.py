from zfel import general_load_bucket

import numpy as np
from scipy import special
import matplotlib.pyplot as plt 
from scipy.signal import find_peaks
import importlib
importlib.reload(general_load_bucket)

def saseTaper(inp_struct):
    '''
    SASE 1D FEL run function
    Input:
    Nruns                       # Number of runs
    npart                       # n-macro-particles per bucket 
    s_steps                     # n-sample points along bunch length
    z_steps                     # n-sample points along undulator
    energy                      # electron energy [MeV]
    eSpread                     # relative rms energy spread [ ]
    emitN                       # normalized transverse emittance [mm-mrad]
    currentMax                  # peak current [Ampere]
    beta                        # mean beta [meter]
    unduPeriod                  # undulator period [meter]
    unduK                       # undulator parameter, K [ ], array of length z_steps
    unduL                       # length of undulator [meter]
    radWavelength               # seed wavelength? [meter], used only in single-freuqency runs
    dEdz                        # rate of relative energy gain or taper [keV/m], optimal~130
    iopt                        # 5=SASE, 4=seeded
    P0                          # small seed input power [W]
    particle_position           # particle information with positions in meter and gamma
    hist_steps                  # number of intervals to generate the histogram of gamma value in a bucket
    constseed                   # whether we want to use constant  random seed for reproducibility, 1 Yes, 0 No
    Output:
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
    gam_out                     # output energy in unit of mc2
    bunching                    # bunching factor
    '''

#eta = (gam - gamavg)/gamavg

    #export variables
    Nruns=inp_struct['Nruns']
    npart=inp_struct['npart']
    s_steps=inp_struct['s_steps']
    z_steps=inp_struct['z_steps']
    energy=inp_struct['energy']
    eSpread=inp_struct['eSpread']
    emitN=inp_struct['emitN']
    currentMax=inp_struct['currentMax']
    beta=inp_struct['beta']
    unduPeriod=inp_struct['unduPeriod']
    unduK=inp_struct['unduK']
    unduL=inp_struct['unduL']
    radWavelength=inp_struct['radWavelength']
    dEdz=inp_struct['dEdz']
    iopt=inp_struct['iopt']
    P0=inp_struct['P0']
    constseed=inp_struct['constseed']
    particle_position=inp_struct['particle_position']
    hist_steps=inp_struct['hist_steps']

    # whether to use constant random seed for reproducibility
    if constseed==1:
        np.random.seed(22)

    # Some constant values
    alfvenCurrent = 17045.0 # Alfven current ~ 17 kA
    mc2 = 0.51099906E6      # Electron rest mass in eV
    c = 2.99792458E8        # light speed in meter/s
    e = 1.60217733E-19      # electron charge in Coulomb
    epsilon_0=8.85418782E-12 #electric constant
    
    gamma0  = energy/mc2                                    # central energy of the beam in unit of mc2
    sigmaX2 = emitN*beta/gamma0                             # rms transverse size, divergence of the electron beam 
    density=currentMax/(e*c*2*np.pi*sigmaX2)                # electron density
    ku=2*np.pi/unduPeriod                                   # undulator wavenumber
    
    #calculating intermediate parameters

    unduJJ  = special.jv(0,np.power(unduK,2)/(4+2*np.power(unduK,2)))\
            -special.jv(1,np.power(unduK,2)/(4+2*np.power(unduK,2)))     # undulator JJ
    
    kappa_1=unduK*unduJJ*e/4/epsilon_0/gamma0               # kappa_1 defined in Kim, Huang, Lindberg
    print(kappa_1)
    Kai= unduK*unduJJ/(2*gamma0**2*mc2)                     # chi_1 defined in Kim, Huang, Lindberg
    print(Kai)

    rho =    np.power(np.power(unduK*unduJJ,2)*np.power(unduPeriod/2/np.pi,2)*\
             (currentMax/alfvenCurrent/(2*sigmaX2)),1/3.0)*(0.5/gamma0)   # Pierce Parameter
    
    resWavelength = unduPeriod*(1+unduK[0]**2/2.0)\
                    /(2.0*gamma0**2)                          # resonant wavelength
    
    Pbeam   = energy*currentMax                           # beam power [W]
    #coopLength = resWavelength/(4*np.pi*rho[0])             # cooperation length
    coopLength = resWavelength/unduPeriod   
    delt  = unduL/z_steps                                   # integration step in z0 ~ 0.1 gain length
    dels  = delt                                            # integration step in s0 must be same as in z0 
    E02   = density*kappa_1[0]*P0*1E-9/Pbeam/Kai[0]         # scaled input power
    
    gbar  = (resWavelength-radWavelength)\
            /(radWavelength)                                # scaled detune parameter
    
    delg  = eSpread                                         # Gaussian energy spread in units of rho 
    
    Ns    = currentMax*unduL/unduPeriod/z_steps\
            *resWavelength/c/e                              # N electrons per s-slice
    
    s = np.arange(1,s_steps+1)*dels*coopLength*1.0e6        # longitundinal steps along beam in meter           
    z = np.arange(1,z_steps+1)*delt                         # longitundinal steps along undulator in meter
    bunchLength=s[-1]*1e-6                                  # beam length in meter
    bunch_steps=np.round(bunchLength/delt/coopLength)       # rms (Gaussian) or half width (flattop) bunch length in s_step
      
     #load buckets
    avg_gamma = np.mean(particle_position[:,1])

    particle_position[:,1] = (particle_position[:,1] - avg_gamma)/(avg_gamma) #normalize gamma

    output_bucket =general_load_bucket.general_load_bucket(npart,Ns, coopLength, particle_position, s_steps, dels, hist_rule='square-root')

    N_real = output_bucket['N_real']
    thet_init = output_bucket['thet_init']
    gam_init = output_bucket['gam_init']


    deta = np.sqrt((1+0.5*unduK[0]**2)/(1+0.5*unduK**2))-1
    print(deta)
    shape = N_real/np.max(N_real)
    
    #shape= 0.5*(np.tanh(10*(np.arange(1,s_steps+1)\
    #       -s_steps/2+bunch_steps)/bunch_steps)\
    #       -np.tanh(10*(np.arange(1,s_steps+1)\
    #       -s_steps/2-bunch_steps)/bunch_steps))            # filling the shape of current and plot it 
    
    plt.figure()
    plt.plot(shape)
   
   
    # sase mode is chosen, go over all slices of the bunch starting from the tail k=1
    if iopt==5: 
        # initialization of variables during the 1D FEL process
        Er=np.zeros((s_steps+1,z_steps+1))
        Ei=np.zeros((s_steps+1,z_steps+1))
        gam=np.zeros((npart,z_steps+1))
        thet_output=np.zeros((npart,z_steps+1))
        thethalf=np.zeros((npart,z_steps+1))
        thet_out=np.zeros((s_steps,1))
        bunching=np.zeros((s_steps,z_steps),dtype=complex)

        for k in range(s_steps):
            Er[k,0] = np.sqrt(E02)                                           # input seed signal
            Ei[k,0] = 0.0
            thet0= thet_init[k,:]
            gam0 = gam_init[k,:]
            gam[:,0] = gam0.T
            thet_output[:,0]=thet0.T
            thethalf[:,0] = thet0.T-2*ku*gam[:,0]*delt/2                            
            thet_out[k,0]=np.mean(thet0.T)
            for j in range(z_steps):                                         # evolve e and gamma in s and t by leap-frog  
                thet = thethalf[:,j]+2*ku*(gam[:,j]+deta[j])*delt/2
                sumsin = np.sum(np.sin(thet))
                sumcos = np.sum(np.cos(thet))
                sinavg = shape[k]*sumsin/npart
                cosavg = shape[k]*sumcos/npart
                Erhalf = Er[k,j]+kappa_1[j]*density * cosavg*dels/2    
                Eihalf = Ei[k,j]-kappa_1[j]*density * sinavg*dels/2   
                thethalf[:,j+1] = thethalf[:,j]+2*ku*(gam[:,j]+deta[j])*delt 
                gam[:,j+1] = gam[:,j]-2*Kai[j]*Erhalf*np.cos(thethalf[:,j+1])*delt\
                             +2*Kai[j]*Eihalf*np.sin(thethalf[:,j+1])*delt  
                thet_output[:,j+1]=thet
                sumsin = np.sum(np.sin(thethalf[:,j+1]))
                sumcos = np.sum(np.cos(thethalf[:,j+1]))
                sinavg = shape[k]*sumsin/npart
                cosavg = shape[k]*sumcos/npart
                Er[k+1,j+1] = Er[k,j]+kappa_1[j]*density *cosavg*dels           # apply slippage condition
                Ei[k+1,j+1] = Ei[k,j]-kappa_1[j]*density *sinavg*dels
                
                bunching[k,j]=np.mean(np.real(np.exp(-1j*thet)))\
                              +np.mean(np.imag(np.exp(-1j*thet)))*1j            #bunching factor calculation
    
        #converting a and gam to field, power and gamavg
        power_s=np.zeros((z_steps,s_steps))
        power_z=np.zeros(z_steps)
        gamavg=np.zeros(z_steps)
        
        for j in range(z_steps):
            for k in range(s_steps):
                power_s[j,k] = (Er[k+1,j]**2+Ei[k+1,j]**2)*(4*np.pi*epsilon_0*sigmaX2*c) 
            power_z[j] = np.sum(Er[:,j]**2+Ei[:,j]**2)*(4*np.pi*epsilon_0*sigmaX2*c)/s_steps 
            
            gamavg[j] = np.sum(gam[:,j+1])/npart                                # average electron energy at every z position
        
        thet_out = thet_output                                                         # don't output phase space
        gam_out = gam
        detune = 2*np.pi/(dels*s_steps)*np.arange(-s_steps/2,s_steps/2+1)
        field = (Er[:,z_steps]+Ei[:,z_steps]*1j)*np.sqrt(Kai[z_steps]/(density*kappa_1[z_steps])*Pbeam)
        field_s = (Er[:,:]+Ei[:,:]*1j)*np.sqrt(Kai/(density*kappa_1)*Pbeam)
        gainLength = unduPeriod/(4*np.pi*rho[0]*np.sqrt(3))                 # gain length from theory
        history={'z':z,'power_z':power_z,'s':s,'power_s':power_s,'field':field,'field_s':field_s,'thet_output':thet_output,'gam':gam,'rho':rho,'detune':detune,'iopt':iopt}
        
    return z,power_z,s,power_s,rho,detune,field,field_s,gainLength,resWavelength,thet_out,gam_out,bunching,history


def plot_log_power_z(history):
    z=history['z']
    power_z=history['power_z']
    plt.figure()
    plt.plot(z,np.log(power_z))
    plt.xlabel('z (m)')
    plt.ylabel('log(P) (W)')

def plot_power_z(history):
    z=history['z']
    power_z=history['power_z']
    plt.figure()
    plt.plot(z,(power_z))
    plt.xlabel('z (m)')
    plt.ylabel('Power (W)')

def plot_power_s(history):
    s=history['s']
    power_s=history['power_s']
    plt.figure()
    for i in range(power_s.shape[0]):
        plt.plot(s,power_s[i,:])
    plt.xlabel('s (um)')
    plt.ylabel('power at different z positions (W)')


def plot_power_s_along_z(history):
    s=history['s']
    power_s=history['power_s']
    z = history['z']
    figs = {}
    axs = {}
    for i in range(power_s.shape[0]):
        if i%4 == 0:
            figs[i] = plt.figure()
            axs[str(z[i])] = figs[i].add_subplot(111)
            axs[str(z[i])].plot(s,power_s[i,:]) #each row is the power at a different z position?
            plt.xlabel('s (m)')
            plt.ylabel('power (W) at z = ' +str(z[i])) #'power at different z positions (W)')
        #plt.show()

def maxP_z(z, power_z):
    dPdz = np.gradient(power_z, z[1]-z[0])
    found = False
    ind = 0

    dPdz_norm = dPdz/power_z[np.argmax(power_z)]
    ind = np.argmax(dPdz_norm)
    '''
    for i in range(dPdz.shape[0]-1):
        if np.sign(dPdz[i+1]) != np.sign(dPdz[i]) and found != True:
            zmax = z[i]
            Pmax = power_z[i]
            ind = i
            found = True
    '''
    #return Pmax, zmax, ind
    return ind

def plot_power_s_at_maxP_z(history):
    s=history['s']
    power_s=history['power_s']
    z = history['z']
    power_z=history['power_z']
    ind = maxP_z(z, power_z)

    plt.figure()
    plt.plot(s,power_s[ind,:])
    plt.xlabel('s (m)')
    plt.ylabel('power (W) at z = ' +str(z[ind]))
    return ind

def delta_s_at_max_P(history):
    s=history['s']
    power_s=history['power_s']
    z = history['z']
    power_z=history['power_z']
    ind = maxP_z(z, power_z)

    found1 = 0
    found2 = 0
    p = power_s[ind,:]
    maxp = p[np.argmax(p)]
    smax = s.shape[0]
    for i in range(smax):
        if p[i] >= .5*maxp and found1 ==0:
            found1 = i
        if p[smax - i - 1] >= 0.5*maxp and found2 == 0:
            found2 = smax-i -1
        if found1 != 0 and found2 !=0:
            break


    return found1, found2

def peaks(history):
    s=history['s']
    power_s=history['power_s']
    z = history['z']
    power_z=history['power_z']
    ind = maxP_z(z, power_z)

    power = power_s[ind,:]
    #peaks at least 1/3 the height of the max power
    p, prop = find_peaks(power, height = power[np.argmax(power)]/3.0)

    return p, prop
    #return # of peaks, and peak position

#will find FWHM of all peaks
#ind index of max z
def peak_width(peaks, ind, history):
    s=history['s']
    power_s=history['power_s']
    width = []
    for i in range(peaks.shape[0]):
        found1 = 0
        found2 = 0
        if i > 0 and s[peaks[i-1]] + width[len(width)-1] > s[peaks[i]]:
            continue #skip if peaks overlap
        for j in range(power_s.shape[0]):
            if peaks[i]-j == 0:
                found1 = 1
                break
            if power_s[ind,peaks[i]-j] <= 0.5*power_s[ind,peaks[i]] and found1 ==0:
                found1 = peaks[i]-j
                break

        for j in range(power_s.shape[0]):
            if peaks[i]+j >= power_s.shape[0] - 1:
                found2 = power_s.shape[0] - 1
                break
            if power_s[ind,peaks[i]+j] <= 0.5*power_s[ind,peaks[i]] and found2 == 0:
                found2 = peaks[i]+j
                break
        print(found2)
        print(found1)
        width.append(s[found2]-s[found1])
    return width





def plot_phase_space(history):
    z=history['z']
    thet_output=history['thet_output']
    gam=history['gam']
    iopt=history['iopt']
    rho=history['rho']
    for j in range(z.shape[0]):
        plt.figure()
        plt.plot(thet_output[:,j],gam[:,j]/rho[j],'.')
        plt.xlabel('theta')
        plt.ylabel('\Delta\gamma/(\gamma rho)')
        '''
        if iopt==4:
            plt.axis([-np.pi,np.pi,-5,5])
        else:
            plt.axis([0,9,-2.5,2.5])
        '''
        plt.title('undulator distance (m) = '+str(z[j]))
       
def plot_pspec(history):
    #need to modify
    field=history['field']
    rho=history['rho']
    detune=history['detune']
    plt.figure()
    fieldFFT = np.fftshift(np.fft(field.T))                  # unconjugate complex transpose by .'
    Pspec=fieldFFT*np.conj(fieldFFT)
    plt.plot(detune*2*rho,Pspec)
    plt.axis([-0.06,+0.01,0,1.2*np.max(Pspec)])
    plt.xlabel('{(\Delta\omega)/\omega_r}')
    plt.ylabel('{output spectral power} (a.u.)')


def plot_norm_power_s(history):
    #need to modify
    power_s=history['power_s']
    z=history['z']
    s=history['s']
    z_steps=z.shape[0]
    s_steps=s.shape[0]
    Y,X = np.meshgrid(z[1:],s);
    p_norm = power_s[1:,:]
    for k in np.arange(0,z_steps-1):
        p_norm[k,:] = p_norm[k,:]/np.max(p_norm[k,:])
    #figure(3)
    #surf(X',Y',p_norm,'EdgeAlpha',0)
    #view(0,90)


def plot_pspec(history):
    #need to modify
    field=history['field']
    rho=history['rho']
    detune=history['detune']
    plt.figure()
    fieldFFT = np.fftshift(np.fft(field.T))                  # unconjugate complex transpose by .'
    Pspec=fieldFFT*np.conj(fieldFFT)
    plt.plot(detune*2*rho,Pspec)
    plt.axis([-0.06,+0.01,0,1.2*np.max(Pspec)])
    plt.xlabel('{(\Delta\omega)/\omega_r}')
    plt.ylabel('{output spectral power} (a.u.)')


#max power = rho*Pbeam

