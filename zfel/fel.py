import numpy as np


def FEL_process_real(npart,
                z_steps,
                kappa_1,
                density,
                Kai,
                ku,
                delt,
                dels,    
                deta,
                thet_init,
                eta_init,
                N_real,
                s_steps,
                E02=0,
                verbose=False):  
    """
    SASE FEL process.
    
    (opt=='sase')
    
    """
    
    if verbose:
        print('FEL_process_real')

    shape = N_real/np.max(N_real) # Scaled current profile shape array. 

    # initialization of variables during the 1D FEL process
    Er  = np.zeros((s_steps+1,z_steps+1))
    Ei  = np.zeros((s_steps+1,z_steps+1))
    eta = np.zeros((npart,z_steps+1))
    thet_output = np.zeros((npart,z_steps+1))
    thethalf    = np.zeros((npart,z_steps+1))
    thet_out    = np.zeros((s_steps,1))
    bunching    = np.zeros((s_steps,z_steps),dtype=complex)
    for k in range(s_steps):
        Er[k,0]  = np.sqrt(E02)                                     # input seed signal
        Ei[k,0]  = 0.0
        thet0    = thet_init[k,:]
        eta0     = eta_init[k,:]
        eta[:,0] = eta0.T
        thet_output[:,0] = thet0.T                                  # eta at j=1
        thethalf[:,0]    = thet0.T-2*ku*eta[:,0]*delt/2             # half back
        thet_out[k,0]    = np.mean(thet0.T)
        for j in range(z_steps):                                    # evolve e and eta in s and t by leap-frog
            thet   = thethalf[:,j]+2*ku*(eta[:,j]+deta[j])*delt/2
            sumsin = np.sum(np.sin(thet))
            sumcos = np.sum(np.cos(thet))
            sinavg = shape[k]*sumsin/npart
            cosavg = shape[k]*sumcos/npart
            Erhalf = Er[k,j]+kappa_1[j]*density * cosavg*dels/2     # minus sign 
            Eihalf = Ei[k,j]-kappa_1[j]*density * sinavg*dels/2               
            thethalf[:,j+1] = thethalf[:,j]+2*ku*(eta[:,j]+deta[j])*delt
            eta[:,j+1] = eta[:,j]-2*Kai[j]*Erhalf*np.cos(thethalf[:,j+1])*delt\
                         +2*Kai[j]*Eihalf*np.sin(thethalf[:,j+1])*delt  #-Eloss*delt  #Eloss*delt to simulate the taper
            
            thet_output[:,j+1] = thet
            sumsin = np.sum(np.sin(thethalf[:,j+1]))
            sumcos = np.sum(np.cos(thethalf[:,j+1]))
            sinavg = shape[k]*sumsin/npart
            cosavg = shape[k]*sumcos/npart
            Er[k+1,j+1] = Er[k,j]+kappa_1[j]*density *cosavg*dels       # apply slippage condition
            Ei[k+1,j+1] = Ei[k,j]-kappa_1[j]*density *sinavg*dels
            bunching[k,j] = np.mean(np.real(np.exp(-1j*thet)))\
                          +np.mean(np.imag(np.exp(-1j*thet)))*1j        #bunching factor calculation

    output = {}
    output['Er'] = Er
    output['Ei'] = Ei  
    output['thet'] = thet_output
    output['eta'] = eta

            
    return output


def FEL_process_complex(npart,
                z_steps,
                kappa_1,
                density,
                Kai,
                ku,
                delt,
                dels,    
                deta,
                thet_init,
                eta_init,
                N_real,
                s_steps,
                E02=0,
                verbose=False):  
    """
    SASE FEL process.
    
    (opt=='sase')
    
    """
    
    if verbose:
        print('FEL_process_complex')

    shape = N_real/np.max(N_real) # Scaled current profile shape array.  

    
    # initialization of variables during the 1D FEL process
    E  = np.zeros((s_steps+1,z_steps+1),dtype=complex) 
    eta = np.zeros((npart,z_steps+1))
    thet_output = np.zeros((npart,z_steps+1))
    thethalf    = np.zeros((npart,z_steps+1))
    #thet_out    = np.zeros((s_steps,1))
    bunching    = np.zeros((s_steps,z_steps),dtype=complex)
    
    # Final bunch
    thet_final  = np.zeros((npart, s_steps))
    eta_final   = np.zeros((npart, s_steps))
    
    # Loop over slices, starting at the back of the bunch
    for k in range(s_steps):
        E[k,0]  = np.sqrt(E02)                                     # input seed signal
        thet0    = thet_init[k,:]
        eta0     = eta_init[k,:]
        eta[:,0] = eta0.T
        thet_output[:,0] = thet0.T                                  # eta at j=1
        thethalf[:,0]    = thet0.T-2*ku*eta[:,0]*delt/2             # half back
        #thet_out[k,0]    = np.mean(thet0.T)

        # Loop over z 
        for j in range(z_steps):                                    # evolve E and eta in s and t by leap-frog

            thet   = thethalf[:,j]+2*ku*(eta[:,j]+deta[j])*delt/2
            
            Ehalf = E[k,j] + kappa_1[j]*shape[k]*density*np.mean(np.exp(-1j*thet))*dels/2  # minus sign, assumes electrons?
                        
            thethalf[:,j+1] = thethalf[:,j]+2*ku*(eta[:,j]+deta[j])*delt  
            
            eta[:,j+1] = eta[:,j]-Kai[j]*delt * 2*np.real(Ehalf*np.exp(1j*thethalf[:,j+1]))
                    
            E[k+1,j+1] = E[k,j]+kappa_1[j]*shape[k]*density*np.mean(np.exp(-1j*thethalf[:,j+1]))*dels  # apply slippage condition
             
            thet_output[:,j+1] = thet # Save for output                
            
            bunching[k,j] = np.mean(np.exp(1j*thet)) #bunching factor calculation    
        
        # Collect final phase space for this slice
        thet_final[:,k] = thet
        eta_final[:,k] = eta[:,-1]
        
        
            
    output = {}
    output['Er'] = np.real(E)
    output['Ei'] = np.imag(E) 
    #output['thet'] = thet_output
    output['thet_final'] = thet_final
    output['eta_final'] = eta_final
    
    # eta history of final slice
    output['theta_final_slice_history'] = thet_output
    output['eta_final_slice_history'] = eta
            
    return output





def final_calc(
               Er,
               Ei,
              # eta,
               s_steps,
               z_steps,
               kappa_1,
               density,
               Kai,
               Pbeam,
               delt,
               dels):
    """
    
    """
 
    #converting a and eta to field, power and etaavg
    power_s = np.zeros((z_steps,s_steps))
    power_z = np.zeros(z_steps)
    #etaavg  = np.zeros(z_steps)
    for j in range(z_steps):
        for k in range(s_steps):
            power_s[j,k] = (Er[k+1,j]**2+Ei[k+1,j]**2)*Kai[j]/(density*kappa_1[j])*Pbeam
        power_z[j] = np.sum(Er[:,j]**2+Ei[:,j]**2)*Kai[j]/(density*kappa_1[j])*Pbeam/s_steps
      #  etaavg[j]  = np.sum(eta[:,j+1])/npart # average electron energy at every z position

    detune   = 2*np.pi/(dels*s_steps)*np.arange(-s_steps/2,s_steps/2+1)
    field    = (Er[:,z_steps]+Ei[:,z_steps]*1j)*np.sqrt(Kai[z_steps-1]/(density*kappa_1[z_steps-1])*Pbeam)
    field_s  = (Er[:,:]+Ei[:,:]*1j)*np.sqrt(np.concatenate((np.array([Kai[0]]),Kai))[np.newaxis,:]/(density*np.concatenate((np.array([kappa_1[0]]),kappa_1))[np.newaxis,:]*Pbeam))
    
    d = {}
    d['power_s'] = power_s
    d['power_z'] = power_z
    
    # Old:
    #pfft = np.fft.fft(field_s[:,1:],axis=0)
    #d['spectrum'] = np.fft.fftshift(np.absolute(pfft)**2)  
    
    # Should apply this over the 0 axis data. 
    def spectrum_from_field(field1):
         return np.abs(np.fft.fftshift(np.fft.fft(field1)))**2
    
    d['spectrum'] = np.apply_along_axis(spectrum_from_field, 0, field_s)
      
    return d
    
