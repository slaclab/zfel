import matplotlib.pyplot as plt
import numpy as np

def plot_log_power_z(history):
    z       = history['z']
    power_z = history['power_z']
    plt.figure()
    plt.yscale('log')
    plt.plot(z, power_z/1e9)
    plt.xlabel('z (m)')
    plt.ylabel('P (GW)')

def plot_power_s(history):
    s = history['s']
    power_s = history['power_s']
    plt.figure()
    for i in range(power_s.shape[0]):
        plt.plot(s,power_s[i,:])
    plt.xlabel('s (um)')
    plt.ylabel('power at different z positions (W)')

def plot_phase_space(history):
    z = history['z']
    thet_output = history['thet_output']
    eta  = history['eta']
    iopt = history['iopt']
    rho  = history['rho']
    for j in range(z.shape[0]):
        plt.figure()
        plt.plot(thet_output[:,j],eta[:,j],'.')
        plt.xlabel('theta')
        plt.ylabel('eta')
        if iopt=='seeded':
            plt.axis([-np.pi,np.pi,-5,5])
        else:
            pass
            #plt.axis([0,9,8400000,8500000])
            #plt.axis([0,9,-2.5,2.5])
        plt.title('undulator distance (m) = '+str(z[j]))
        #pause(.02)

def plot_pspec(history):
    #need to modify
    field  = history['field']
    rho    = history['rho']
    detune = history['detune']
    plt.figure()
    fieldFFT = np.fftshift(np.fft(field.T))                  # unconjugate complex transpose by .'
    Pspec=fieldFFT*np.conj(fieldFFT)
    plt.plot(detune*2*rho,Pspec)
    plt.axis([-0.06,+0.01,0,1.2*np.max(Pspec)])
    plt.xlabel('{(\Delta\omega)/\omega_r}')
    plt.ylabel('{output spectral power} (a.u.)')


def plot_norm_power_s(history):
    #need to modify
    power_s = history['power_s']
    z = history['z']
    s = history['s']
    z_steps = z.shape[0]
    s_steps = s.shape[0]
    Y,X = np.meshgrid(z[1:],s);
    p_norm = power_s[1:,:]
    for k in np.arange(0,z_steps-1):
        p_norm[k,:] = p_norm[k,:]/np.max(p_norm[k,:])


def plot_current(history):
    shape = history['shape']
    plt.figure()
    plt.plot(shape)
    plt.title('Current')

