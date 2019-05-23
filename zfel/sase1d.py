import numpy as np
import scipy
import load_bucket
import matplotlib.pyplot as plt 

def sase(inp_struct):
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
    unduK                       # undulator parameter, K [ ]
    unduL                       # length of undulator [meter]
    radWavelength               # seed wavelength? [meter], used only in single-freuqency runs
    dEdz                        # rate of relative energy gain or taper [keV/m], optimal~130
    iopt                        # 5=SASE, 4=seeded
    P0                          # small seed input power [W]
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

    # whether to use constant random seed for reproducibility
    if constseed==1:
        np.random.seed(22)

    # Some constant values
    alfvenCurrent = 17045.0 # Alfven current ~ 17 kA
    mc2 = 510.99906E-3      # Electron rest mass in MeV
    c = 2.99792458E8        # light speed in meter
    e = 1.60217733E-19      # electron charge in Coulomb

    #calculating intermediate parameters
    unduJJ  = scipy.special.jv(0,unduK**2/(4+2*unduK**2))\
              -scipy.special.jv(1,unduK**2/(4+2*unduK**2))  # undulator JJ
    gamma0  = energy/mc2                                    # central energy of the beam in unit of mc2
    sigmaX2 = emitN*beta/gamma0                             # rms transverse size, divergence of the electron beam
    rho     = (0.5/gamma0)*((currentMax/alfvenCurrent)\
              *(unduPeriod*unduK*unduJJ/(2*np.pi))**2\
              /(2*sigmaX2))**(1/3)                          # FEL Pierce parameter
    resWavelength = unduPeriod*(1+unduK**2/2.0)\
                    /(2*gamma0**2)                          # resonant wavelength
    rhoPbeam   = rho*energy*currentMax/1000.0               # rho times beam power [GW]
    coopLength = resWavelength/(4*np.pi*rho)                # cooperation length
    gainLength = unduPeriod/(4*np.pi*rho)                   # rough gain length
    #cs0  = bunchLength/coopLength                           # bunch length in units of cooperation length     
    z0    = unduL/gainLength                                # wiggler length in units of gain length
    delt  = z0/z_steps                                      # integration step in z0 ~ 0.1 gain length
    dels  = delt                                            # integration step in s0 must be same as in z0 
    a02   = P0*1E-9/rhoPbeam                                # scaled input power
    gbar  = (resWavelength-radWavelength)\
            /(radWavelength*rho)                            # scaled detune parameter
    delg  = eSpread/rho                                     # Gaussian energy spread in units of rho 
    Ns    = currentMax*unduL/unduPeriod/z_steps\
            *resWavelength/c/e                              # N electrons per s-slice [ ]
    Eloss = -dEdz*1E-3/energy/rho*gainLength                # convert dEdz to alpha parameter
    s = np.arange(1,s_steps+1)*dels*coopLength*1.0e6        # longitundinal steps along beam in meter           
    z = np.arange(1,z_steps+1)*delt*gainLength              # longitundinal steps along undulator in meter

    bunchLength=s[-1]*1e-6/2.5                              # beam length in meter
    bunch_steps=np.round(bunchLength/delt/coopLength)       # rms (Gaussian) or half width (flattop) bunch length in s_step
    shape = np.zeros((1,s_steps))                           # initialization of the beam current shape
    shape= 0.5*(np.tanh(10*(np.arange(1,s_steps+1)\
           -s_steps/2+bunch_steps)/bunch_steps)\
           -np.tanh(10*(np.arange(1,s_steps+1)\
           -s_steps/2-bunch_steps)/bunch_steps))            # filling the shape of current and plot it
    plt.plot(shape)

    # sase mode is chosen, go over all slices of the bunch starting from the tail k=1
    if iopt==5: 
        # initialization of variables during the 1D FEL process
        ar=np.zeros((s_steps+1,z_steps+1))
        ai=np.zeros((s_steps+1,z_steps+1))
        gam=np.zeros((npart,z_steps+1))
        thet_output=np.zeros((npart,z_steps+1))
        thethalf=np.zeros((npart,z_steps+1))
        thet_out=np.zeros((s_steps,1))
        bunching=np.zeros((s_steps,z_steps),dtype=complex)
        for k in range(s_steps):
            ar[k,0] = np.sqrt(a02)                                              # input seed signal
            ai[k,0] = 0.0
            [thet0,gam0] = load_bucket.load_bucket(npart,gbar,delg,iopt,Ns)     # load each bucket
            gam[:,0] = gam0.T
            thet_output[:,0]=thet0.T                                                   # gamma at j=1
            thethalf[:,0] = thet0.T-gam[:,0]*delt/2                             # half back
            thet_out[k,0]=np.mean(thet0.T)
            for j in range(z_steps):                                            # evolve e and gamma in s and t by leap-frog
                thet = thethalf[:,j]+gam[:,j]*delt/2
                sumsin = np.sum(np.sin(thet))
                sumcos = np.sum(np.cos(thet))
                sinavg = shape[k]*sumsin/npart
                cosavg = shape[k]*sumcos/npart
                arhalf = ar[k,j]+cosavg*dels/2
                aihalf = ai[k,j]-sinavg*dels/2
                thethalf[:,j+1] = thethalf[:,j]+gam[:,j]*delt
                gam[:,j+1] = gam[:,j]-2*arhalf*np.cos(thethalf[:,j+1])*delt\
                             +2*aihalf*np.sin(thethalf[:,j+1])*delt-Eloss*delt  #Eloss*delt to simulate the taper
                thet_output[:,j+1]=thet
                sumsin = np.sum(np.sin(thethalf[:,j+1]))
                sumcos = np.sum(np.cos(thethalf[:,j+1]))
                sinavg = shape[k]*sumsin/npart
                cosavg = shape[k]*sumcos/npart
                ar[k+1,j+1] = ar[k,j]+cosavg*dels                               # apply slippage condition
                ai[k+1,j+1] = ai[k,j]-sinavg*dels
                bunching[k,j]=np.mean(np.real(np.exp(-1j*thet)))\
                              +np.mean(np.imag(np.exp(-1j*thet)))*1j            #bunching factor calculation
    
        #converting a and gam to field, power and gamavg
        power_s=np.zeros((z_steps,s_steps))
        power_z=np.zeros(z_steps)
        gamavg=np.zeros(z_steps)
        for j in range(z_steps):
            for k in range(s_steps):
                power_s[j,k] = (ar[k+1,j]**2+ai[k+1,j]**2)*rhoPbeam
            power_z[j] = np.sum(ar[:,j]**2+ai[:,j]**2)/s_steps*rhoPbeam
            gamavg[j] = np.sum(gam[:,j+1])/npart                                # average electron energy at every z position
            thet_out=0                                                          # don't output phase space
            gam_out=0
        detune = 2*np.pi/(dels*s_steps)*np.arange(-s_steps/2,s_steps/2+1)
        field = (ar[:,z_steps]+ai[:,z_steps]*1j)*np.sqrt(rhoPbeam)
        field_s = (ar[:,:]+ai[:,:]*1j)*np.sqrt(rhoPbeam)

        history={'z':z,'power_z':power_z,'s':s,'power_s':power_s,'field':field,'field_s':field_s,'thet_output':thet_output,'gam':gam,'rho':rho,'detune':detune,'iopt':iopt}




    else:            # seeded mode then
        thet_out=np.zeros((npart,z_steps))
        gam_out=np.zeros((npart,z_steps))
        power_s=np.zeros((z_steps,1))
        power_z=np.zeros(z_steps)
        gamavg=np.zeros(z_steps)



        delth=delt/2.0
        [thet,gam] = load_bucket.load_bucket(npart,gbar,delg,iopt,Ns)
        ar=np.sqrt(a02)      # initial seed signal
        print(ar)
        ai=0.0
     
        for j in range(z_steps):
        #            first RK step 
            thet1=thet+gam*delth
            gam1=gam-2*ar*np.cos(thet)*delth+2*ai*np.sin(thet)*delth-Eloss*delth
            sinavg=np.sum(np.sin(thet))/npart
            cosavg=np.sum(np.cos(thet))/npart
            ar1=ar+cosavg*delth
            ai1=ai-sinavg*delth
        #            second RK step
            thet2=thet+gam1*delth
            gam2=gam-2*ar1*np.cos(thet1)*delth+2*ai1*np.sin(thet1)*delth-Eloss*delth
            sinavg=np.sum(np.sin(thet1))/npart
            cosavg=np.sum(np.cos(thet1))/npart
            ar2=ar+cosavg*delth
            ai2=ai-sinavg*delth
        #            third RK step
            thet3=thet+gam2*delt
            gam3=gam-2*ar2*np.cos(thet2)*delt+2*ai2*np.sin(thet2)*delt-Eloss*delt
            sinavg=np.sum(np.sin(thet2))/npart
            cosavg=np.sum(np.cos(thet2))/npart
            ar3=ar+cosavg*delt
            ai3=ai-sinavg*delt
        #            fourth RK stfigure(3)
            thet4=thet+gam3*delth
            gam4=gam-2*ar3*np.cos(thet3)*delth+2*ai3*np.sin(thet3)*delth-Eloss*delth
            sinavg=np.sum(np.sin(thet3))/npart
            cosavg=np.sum(np.cos(thet3))/npart
            ar4=ar+cosavg*delth
            ai4=ai-sinavg*delth
        #            add them up
            thet=thet1/3+thet2*2/3+thet3/3+thet4/3-thet*2/3;   
            gam=gam1/3+gam2*2/3+gam3/3+gam4/3-gam*2/3
            ar=ar1/3+ar2*2/3+ar3/3+ar4/3-ar*2/3   #-2/3*ar because already added 5/3*ar
            ai=ai1/3+ai2*2/3+ai3/3+ai4/3-ai*2/3
        # thet_out and gam_out can be use for longitudinal phase space display 
        # to output phase between -pi to +pi
            for i in range(npart):
                while thet[i]>np.pi:
                    thet[i]=thet[i]-2*np.pi
                while thet[i]<-np.pi:
                    thet[i]=thet[i]+2*np.pi
            thet_out[:,j]=thet.T
            gam_out[:,j]=gam.T
        #           output
            power_z[j]=(ar**2+ai**2)*rhoPbeam
            power_s[j,:] = power_z[j]
            gamavg[j] = np.sum(gam)/npart         # beam power loss
            detune = 0
            field = 0
            field_s = 0
            bunching= 0
        history={'z':z,'power_z':power_z,'s':s,'power_s':power_s,'field':field,'field_s':field_s,'thet_output':thet_out,'gam':gam_out,'rho':rho,'detune':detune,'iopt':iopt}

    return z,power_z,s,power_s,rho,detune,field,field_s,gainLength,resWavelength,thet_out,gam_out,bunching,history


def plot_log_power_z(history):
    z=history['z']
    power_z=history['power_z']
    plt.figure()
    plt.plot(z,np.log(power_z*1E9))
    plt.xlabel('z (m)')
    plt.ylabel('log(P) (W)')

def plot_power_s(history):
    s=history['s']
    power_s=history['power_s']
    plt.figure()
    for i in range(power_s.shape[0]):
        plt.plot(s,power_s[i,:])
    plt.xlabel('s (m)')
    plt.ylabel('power at different z positions (W)')

def plot_phase_space(history):
    z=history['z']
    thet_output=history['thet_output']
    gam=history['gam']
    iopt=history['iopt']
    for j in range(z.shape[0]):
        plt.figure()
        plt.plot(thet_output[:,j],gam[:,j],'.')
        plt.xlabel('theta')
        plt.ylabel('\Delta\gamma/(\gamma rho)')
        if iopt==4:
            plt.axis([-np.pi,np.pi,-5,5])
        else:
            plt.axis([0,9,-2.5,2.5])
        plt.title('undulator distance (m) = '+str(z[j]))
        #pause(.02)

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







