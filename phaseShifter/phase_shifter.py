# Copyright (C) 2017  Gordon, Sam <sbg2133@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.signal import spectrogram

################################################
# W-Band NbTiN Microstrip phase shifter circuit
# 75 - 110 GHz, 2.4 - 7 mm
# All units are SI, unless otherwise noted
################################################

# To minimize vPhase, want d, t << lamda
# Note: "Because microstrip line is not a true TEM line, its propagation
# constant (B = sqrt(er) * k_0  is not a linear function of frequency, meaning
# that the effective dielectric constant varies with frequency." Pozar pg. 150

# kinetic inductance units: [Ohm * s], or [kg / m^4 * e^2 ]

###########################################################
# Physical constants
boltz = 1.38062e-23 # boltzmann constant, k [J/K]
hBar = 1.054e-34 # reduced planck constant, hbar [J*s / rad]
h = hBar * 2.*np.pi # planck constant, h [J*s]
e0 = 8.854e-12 # permitivity of free space, [Farad/m]
mu0 = 4.0*np.pi * 1.0e-7 # permeability of free space, [H/m]
c = 3.0e8 # speed of light, [m/s]
#############################################################
# Geometric parameters
l_max = 503.0e-3 # length of longest delay line, [m]
d = 20.0e-6 # dielectric thickness, [m]
t = 20.0e-9 # approx. thickness of superconductor, [m]
w = 3.0e-6 # approx. thickness of Becky's substrate, [m]
tbecky = 40.0e-9
###############################################################
# Material parameters
Tc = 12.0 # approx. critical temperature for NbTiN, [K]
rho = 100.0e-8 # normal state resistivity of NbTiN, [Ohm * m]
er = 1.0 # Assumed relative permitivity
lonDep0 = 200.0e-9 # London penetration depth, NbTiN [m]
gap = (1.76 * boltz * Tc) # Superconducting energy gap, [J]
dnu_nu = 2.5e-4 # Fractional frequency shift, Becky's resonator
#############################################################
# Current parameters
Istar0 = 0.4e-3 # 0.4 mA, critical current in Becky's resonator w/ w ~ 10 um
jStar0 = Istar0 / ( t * 10.0e-6) # critical current density, [A/m^2]
Istar = jStar0 * w * t # Critical current, amps
#############################################################
# Measurement variables
Tsys = 4.0 # system temperature, [K]
pumpFreq = 11.56e9 # Pump frequency, Nature paper, [Hz]
freqs = np.arange(75., 112.5, 2.5)*1.0e9
I = np.arange(1.0e-10, Istar, 1.0e-8)
#I = np.arange(Istar/1000.,np.sqrt(4)* Istar, Istar/1000.)
##############################################################

def lonDepth(T):
    """London penetration depth, in units of m
       @param T: Temperature of circuit, in K"""
    return lonDep0 * ( 1.0  - ((T/Tc)**4)**(-1/2))

def nu_0(L, C):
    """Resonant frequency
       @param L: Inductance
       @param C: Capacitance"""
    return 1. / (2.*np.pi * np.sqrt(L * C))

#def dnu_nu( L_tot, L_k0, nu_0 ):
#    """Phil's slides, maximum frequency shift before loss"""
#    return( alpha / 2.) * (4./9.)*( I_max**2. / Istar**2.)

def e_nu(freqs, e_e, Z, d ):
    """Dielectric constant as a function of frequency (Pozar 3.200)
        @param freqs: array of frequencies, in GHz
        @param e_e: effective dielectric constant
        @param Z: Impedance
        @param d: dielectric thickness"""
    g = 0.6 + (0.009 * Z)
    nu_p = Z / ( 8. * np.pi * (d / 1.0e-2))
    G_nu = g * ( freqs / nu_p )**2.
    e_nu = er - ( (e_r - e_e ) / (1. + G_nu))
    return e_nu, G_nu

def plot_enu():
    nu = np.arange(0., 20.)
    Z_0 = 25.
    ee = e_e(w,d)
    enu, gnu = e_nu(freqs, ee,  Z_0, d)
    print 'Z =',  Z_0
    print 'er=', e_r
    print 'ee0 =', ee
    return freqs, enu, gnu

def e_e(w, d):
    """Effective dielectric constant
       See Pozar, section 3.8: 1 < e_e < e_r
       @param w: width of microstrip
       @param d: dielectric thickness"""
    return ( (er + 1.0) / 2.0 ) + (( (e_r - 1.0) / 2.0) /\
                     (np.sqrt( 1.0 + 12.0 * (d / w))))

def beta(e_e, lam):
    """Beta, the propagation constant
       @param e_e: effective dielectric constant
       @param lam: wavelength"""
    return (2.0*np.pi / lam) * np.sqrt(e_e)

########################################
# Functions that compute inductances
#########################################

def LkRho():
    """Intrinsic kinetic inductance, per unit length"""
    return (hBar * rho) / (np.pi * gap * w * t)

### From Simon Doyle's thesis ###
def LkThin():
    """Intrinsic KI for thin films, w << lam << t """
    x = (t / (2.* lonDepth0))
    return ( (mu0 * lonDepth0) / (4. * w)) * ( (1./np.tanh(x)) +\
              (x * np.square(1./np.sinh(x))))

### Simon Doyle's thesis ###
def LGeo_SD():
    """Geometric inductance, per unit length for thin film"""
    x = (t / (2. * lonDepth0))
    return ( (mu0 * lonDepth0) / (4. * w)) * ( (1./np.tanh(x)) -\
              (x * np.square(1./np.sinh(x))))

def plot_L():
    #t = np.arange(0.0, 50.0)*1.0e-9 # nm
    tot = LkThin() + LGeo()
    #tot = ( (mu0 * lonDepth0) / (2.)) * (1./np.tanh(x))
    plt.plot(t*1.0e9, LkThin(), label = r'L${_k}$')
    plt.plot(t*1.0e9, LGeo(), label = r'L${_m}$')
    plt.xlabel('Film thickness [nm]')
    plt.ylabel('Inductance')
    plt.title('Fractions for fixed w')
    plt.legend(loc = 'upper right')
    return

def L_m():
    """Geometric Inductance per unit length"""
    return mu0 * (d / w )

def LGeo(t):
    """Geometric inductance per unit length, from Harshad's Sonnet sim"""
    return 1.11e-6 * t / (20.0e-9)

def LkLin():
    """Intrinsic KI per unit length, fn of lonDepth
       From Pond paper, assume KI >> L_geo"""
    return (mu0 * lonDep0**2.) / (w * t)  # * ((1. - (T/Tc)**4.)**-1. )

def LkNonlin(I_in):
    """Nonlinear KI, fn of LkRho
       @param I_in: Current, Amps"""
    return LkRho()*( 1. + ( I_in**2. / Istar**2. ))

def Ltot(I_in):
    """Total inductance, KI + Lgeo
       @param I_in: Current, Amps"""
    Lk = LkNonlin(I_in)
    Lm = LGeo(t)
    return ( Lk + Lm ) * l_max

def LkFrac(I_in):
    """Alpha, the inductance fraction
       @param I_in: Current, Amps"""
    tot = LkRho() + LGeo(t)
    return LkNonlin(I_in) / tot

def plot_Lk():
    tot = Ltot(I_in)
    plt.plot(I_in*1.0e3, tot)
    plt.plot(I_in*1.0e3, tot)
    plt.title('Kinetic Inductance fraction vs. DC Bias Current')
    plt.xlabel('Current (uA)')
    plt.ylabel('Lk fraction')
    plt.show()
    return

############################
# Capacitance
############################

def Cap_unused():
    """Capacitance per unit length"""
    return (e0 * er * w )/ d

def Cap():
    """Capacitance per unit length
       Determined from two Sonnet simulations: Microstrip with KI/square,
           and PEC (no KI).
       For t = 20 nm"""
    return 8.38e-11

###################
# Phase velocities
###################

def v_0():
    """From Pond, assume, KI >> Lgeo"""
    return (1./lonDep0) * np.sqrt( (d * t )/ (2. * mu0 * e0 * er))

def vPhase(L, C):
    """Phase velocity,
       @param L: inductance
       @param C: capacitance"""
    return 1.0 / np.sqrt( L * C  )

####################
# Impedance
####################

def Z(I_in):
    """Impedance
       @param I_in: current in amps"""
    Z =  np.sqrt((Ltot(I_in)/l) / (Cap()))
    plt.title('Strip Impedence vs. Current')
    plt.xlabel('Current (mA)')
    plt.ylabel('Z (I) (Ohms)')
    plt.plot(I_in*1.0e3, Z)
    plt.show()
    return Z

###########################
# Phase (phi) and phase shift
############################

def phi(freqs,l, v):
    """Phase shift
       @param freqs: array of freqs [Hz]
       @param l: propagation length, [m]
       @param v: phase velocity, [m/s]"""
    phi = (2. * np.pi * freqs * l) / v
    return phi

def pShift(freqs, I_in, plot = False):
    """Phase shift as a function of current for one or more frequencies.
       @param freqs: array of freqs [Hz]
       @param I_in: An array of currents"""
    if plot:
        plt.figure(figsize=(10.24, 7.68), dpi = 100)
        ax = plt.gca()
    v0 = vPhase(Ltot(0)/l_max, Cap())
    v = vPhase(Ltot(I_in)/l_max, Cap())
    pShift = np.zeros((len(freqs), len(I_in)))
    phi0 = phi(freqs, l_max, v0)
    for i in range(len(freqs)):
        phis = phi(freqs[i], l_max, v)
        pShift[i] = phis - phi0[i]
    if plot:
        colormap = plt.cm.OrRd
        plt.gca().set_color_cycle([colormap(i) for\
                        i in np.linspace(0.3, 1.0, len(freqs))])
        ax.tick_params(axis='x', labelsize=18)
        ax.tick_params(axis='y', labelsize=18)
        [i.set_linewidth(3.) for i in ax.spines.itervalues()]
        for i in range(len(freqs)):
            phis = phi(freqs[i], l_max, v)
            pShift[i] = phis - phi0[i]
            lab = str(freqs[i]*1.0e-9)
            ax.plot(I_in*1.0e3, -pShift[i], label = lab)
    if plot:
        ax.set_title('Phase Shifter: Predicted current-dependent shift',\
         size = 18, fontweight='bold')
        ax.set_xlabel('Current (mA)', fontsize = 18, fontweight='bold')
        ax.set_ylabel(r'$\Delta$$\phi$ (rad)', fontsize = 18,\
         fontweight='bold')
        ax.legend(loc = 'lower left', fontsize = 18, ncol = 2)
        plt.grid()
        plt.savefig('./pShift.png', bbox_inches='tight')
    return pShift

############################
# Functions that compute gain
#############################

def quadGain(freqs, I_in):
    """Quadratic gain, [dB], for all currents in I array (global)
       @param phi: phase, rad
       @param dnu_nu: fractional frequency shift"""
    dB = np.zeros((len(freqs), len(I_in)))
    p_shifts = pShift(freqs, I_in, plot = False)
    for i in range(len(freqs)):
        G = 1. + (dnu_nu * p_shifts[i])**2.
        dB[i] = 10.*np.log10(G)
    return dB

def quadGainOneCurrent(freqs, I_in):
    """Quadratic gain [dB], for a single current"""
    p_shift = pShift(freqs, I_in, plot = False)
    for i in range(len(freqs)):
        G = 1. + (dnu_nu * p_shifts)**2.
    return 10.*np.log10(G)

def plotPshiftGain(freqs, I_in):
    """Plots the current dependent phase shift and quadratic gain for
       each frequency in freqs"""
    plt.figure(figsize=(10.24, 7.68), dpi = 100)
    ax = plt.gca()
    shift = pShift(freqs, I_in, plot = False)
    gain = quadGain(freqs, I_in)
    colormap = plt.cm.OrRd
    plt.gca().set_color_cycle([colormap(i)\
                        for i in np.linspace(0.4, 1.0, len(freqs))])
    plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
    [i.set_linewidth(3.) for i in ax.spines.itervalues()]
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    [i.set_linewidth(3.) for i in ax.spines.itervalues()]
    [ax.plot(I*1.0e3, -shift[i], linewidth = 2,\
                   label = str(freqs[i]*1.0e-9)) for i in range(len(freqs))]
    ax2 = ax.twinx()
    ax2.set_ylabel('Quadratic Gain [dB]', fontsize = 18, fontweight='bold')
    [i.set_linewidth(3.) for i in ax2.spines.itervalues()]
    colormap = plt.cm.Blues
    plt.gca().set_color_cycle([colormap(i)\
                          for i in np.linspace(0.4, 1.0, len(freqs))])
    plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax2.tick_params(axis='x', labelsize=18)
    ax2.tick_params(axis='y', labelsize=18)
    [ax2.plot(I*1.0e3, gain[i], linewidth = 2) for i in range(len(freqs))]
    ax.set_title('Phase Shift and Quadratic Gain v. bias current',\
          size = 18, fontweight='bold')
    ax.legend(loc ='lower left', fontsize = 18, ncol = 2)
    plt.grid()
    plt.savefig('./pShiftGain.png',bbox_inches='tight')
    return

######################################
# dielectric loss, noise temperatures
######################################

def alpha_d_poz(freqs, e_e):
    """Dielectric loss coefficient, from Pozar
       tan_del = e" / e' taken to be 2.0e-4, from Pond paper
       Here er ( e_e - 1) / (sqrt(e_e) (e_r - 1) is a 'filling factor'
       to account for the wave being partially in air and partially in
          the dielectric
       @param freqs: array of frequencies, [Hz]
       @param e_e: effective permitivity"""
    return ((2.0 * np.pi * freqs/ c) * er * (e_e - 1.0) * tan_del) /\
         (2.0 * np.sqrt(e_e) * (e_r - 1.0)) # Np / m

def Tdiel(freqs, v, Tsys, l):
    """Noise temperature of dielectric
       @param freqs: array of freqs, [Hz]
       @param v: phase velocity
       @param Tsys: System temp
       @param l: propagation length"""
    alpha_d = (((2.*np.pi * freqs)/ (2 * v)) * ( 1.0e-5) )   # dB/m
    loss = 1.0 - np.exp(-2.0* alpha_d * l)
    return loss * Tsys

def Tquant(freqs):
    """Quantum noise temperature
       @param freqs: array of freqs [Hz]"""
    return (h * freqs) / (2. *boltz)

def noiseComp(I_in):
    Tdie = np.zeros((len(freqs), len(I_in)))
    Tq = np.zeros((len(freqs),len(I_in)))
    v = vPhase(Ltot(I_in)/l, Cap())
    p = I**2. * Z(I_in)
    for i in range(len(freqs)):
        Tq[i] = Tquant(freqs[i])
        Tdie[i] = Tdiel(freqs[i], v, Tsys, l)
        lab =( str( freqs[i]*1.0e-9 ))
        plt.plot(p * 1.0e6, Tq[i]/Tdie[i], label = lab)
    plt.title(r'N$_{Q}$ / N$_{d}$, T$_{sys}$ = 4 K')
    plt.ylabel(r'N$_{Q}$ / N$_{d}$')
    plt.xlabel( 'Pump Power (uW)')
    plt.legend(loc = 'upper right')
    plt.show()
    return

def maxBias(del_T):
    """DC current handling capacity, delta T < 100 K,
       K Silicon Dioxide taken as 1.4 W / m * K"""
    return w * np.sqrt(  149. * t * del_T / ( rho * d) )

def plot_Z_w():
    t = np.arange(0.0, 51.0)*1.0e-9
    w = np.array([2, 4, 8, 10, 12])*1.0e-6
    Zs = np.zeros((len(w), len(t)))
    for i in range(len(w)):
        Zs[i] = Z(w[i], t, d)
    vPhase = v_0(t, d)
    fig, ax1 = plt.subplots()
    plt.title(r'Z$_0$ vs. v$_p$')
    ax1.tick_params(axis='x', pad=15)
    ax1.set_xlabel( 'Film Thickness (nm)')
    ax1.set_ylabel(r'Z$_0$ (ohms)', color='b')
    ax1.yaxis.set_ticks(np.arange(0.0, 501.0, 25.0))
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Normalized Phase Velocity (v / c)', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    for i in range(len(w)):
        lab = 'w = '  + str( w[i] * 1.0e6 ) + ' um'
        ax1.plot(t*1.0e9, Zs[i], label = lab)
        #ax1.text( 2.0, Zs[i][0], lab)
    ax2.plot(t*1.0e9, vPhase / c, 'r--')
    ax1.legend(loc = 'upper middle')
    plt.show()
    return

#######################################
# Four Wave Mixing --> Exponential Gain
#######################################

def waveNum(w, I_in, l):
    """Wavenumber
       @param w: microstrip width
       @param I_in: current
       @param l: propagation length"""
    k = w / vPhase(Ltot(I_in)/l, Cap())
    return k

def parampGain(I_in, l, sigFreq, pumpFreq, plot = True):
    """Computes the exponential gain as a function of propagation length
       along the length of the microstrip line, for a single signal and
       pump frequency.
       @param I_in: current
       @param l: length of line
       @param sigFreq: signal frequency
       @param pumpFreq: pump frequency
       @param plot: If True, plot results"""
    if plot:
        plt.figure(figsize = (10.24, 7.68), dpi = 100)
        ax = plt.gca()
    wp = 2.*np.pi * pumpFreq
    ws = 2.*np.pi * sigFreq
    wi = (2*wp) - ws
    alpha = LkFrac(I_in)
    Cp =  1.j * alpha * waveNum(wp, I_in, l) / (2.*Istar**2)
    Cs =  1.j * alpha * waveNum(ws, I_in, l) / (2.*Istar**2)
    Ci =  1.j * alpha * waveNum(wi, I_in, l) / (2.*Istar**2)
    dz = np.arange(0.001, l, 1.0e-4)
    Ap = np.zeros(len(dz), dtype = 'complex')
    As = np.zeros(len(dz), dtype = 'complex')
    Ai = np.zeros(len(dz), dtype = 'complex')
    dAp = 1.0e-3*I_in
    dAs = 1.0e-7*I_in
    dAi = 1.0e-10*I_in
    for i in range(len(dz)):
        dAp += Cp * dz[i] * (( np.abs(dAp)**2 + 2.*np.abs(dAs)**2 +\
               2.*np.abs(dAi)**2 ) * dAp + 2.* dAs * dAi * np.conjugate(dAp)*\
               np.exp((1.j * -waveNum(wp, I_in, l) * alpha * np.abs(dAp)**2 *\
               dz[i]) / (2.*Istar**2)) )
        dAs += Cs * dz[i] * (( np.abs(dAs)**2 + 2.*np.abs(dAi)**2 +\
               2.*np.abs(dAp)**2 ) * dAs +  np.conjugate(dAi) * dAp**2 *\
               np.exp((1.j * -waveNum(wp, I_in, l) * alpha * np.abs(dAp)**2 *\
               dz[i]) / (2.*Istar**2)))
        dAi += Ci * dz[i] * (( np.abs(dAi)**2 + 2.*np.abs(dAs)**2 +\
               2.*np.abs(dAp)**2 ) * dAi +  np.conjugate(dAs) * dAp**2 *\
               np.exp((1.j * -waveNum(wp, I_in, l) * alpha * np.abs(dAp)**2 *\
               dz[i]) / (2.*Istar**2)))
        Ap[i] = dAp
        As[i] = dAs
        Ai[i] = dAi
    Pp = (Ap**2) * 50.
    Sp = (As**2) * 50.
    Ip = (Ai**2) * 50.
    Ap_dBm = 10.*np.log10((Pp*1.0e3))
    As_dBm = 10.*np.log10((Sp*1.0e3))
    Ai_dBm = 10.*np.log10((Ip*1.0e3))
    Gain = np.abs(np.abs(As_dBm[-10]) - np.abs(As_dBm[10]))
    if plot:
        ax.tick_params(axis='x', labelsize=18)
        ax.tick_params(axis='y', labelsize=18)
        [i.set_linewidth(3.) for i in ax.spines.itervalues()]
        ax.plot(dz*1.0e3, Ap_dBm, color = 'b', label = 'Pump', linewidth = 3)
        ax.plot(dz*1.0e3, As_dBm, color = 'r', label = 'Signal', linewidth = 3)
        ax.plot(dz*1.0e3, Ai_dBm, color = 'c', label = 'Idler', linewidth = 3)
        plt.ylim((-250, -80))
        plt.ylabel('Power (dBm)', fontsize = 18, fontweight='bold')
        plt.xlabel('Propagation length (mm)', fontsize = 18,\
                              fontweight='bold')
        plt.title(r'Signal Gain = ' + str(np.round(Gain, 2)) +\
                      'dB, $\omega_{p}$ = ' + str(pumpFreq/1.0e9)\
                      + ', GHz, $\omega_{s}$ = ' + str(sigFreq/1.0e9) +\
                      ' GHz', fontsize = 18, fontweight='bold')
        plt.legend(loc = 'lower right', fontsize = 18)
        plt.grid()
        plt.savefig('./parametricGain.png', bbox_inches='tight')
    return Gain

def parampFreqResponse(l, pumpFreq):
    """Plots the frequency response, in dB, of the paramp
       @param l: line length, m
       @param pumpFreq: pump frequency, [GHz]"""
    plt.figure(figsize=(10.24, 7.68), dpi = 100)
    ax = plt.gca()
    freqs = np.linspace(70., 115., 150)*1.0e9
    gain = np.zeros(len(freqs))
    frac_freq = np.zeros(len(freqs))
    for i in range(len(freqs)):
    	gain[i] = parampGain(Istar, l, freqs[i],\
                     pumpFreq, plot = False)
    	frac_freq[i] = freqs[i] / pumpFreq
    ax.axvline(pumpFreq/1.0e9, 0, 1, color = 'r', linestyle = '--')
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    [i.set_linewidth(3.) for i in ax.spines.itervalues()]
    ax.plot(np.round(freqs/1.0e9, 2), gain, linewidth = 3)
    plt.title(r'W-Band Paramp Freq Response, f$_{pump}$ = ' + str(pumpFreq/1.0e9)\
                           + ' GHz', fontsize = 18, fontweight = 'bold')
    plt.ylabel('Gain (dB)', fontsize = 18, fontweight = 'bold')
    plt.xlabel('frequency (GHz)', fontsize = 18, fontweight = 'bold')
    plt.grid()
    plt.savefig('./freqResponse.png', bbox_inches='tight')
    return

def pShifterFreqResponse(l, I, pumpFreq):
    """Plots the frequency response, in dB, of the paramp
       @param l: line length, m
       @param pumpFreq: pump frequency, [GHz]"""
    plt.figure(figsize=(10.24, 7.68), dpi = 100)
    ax = plt.gca()
    freqs = np.linspace(70., 115., 150)*1.0e9
    gain = quadGain(freqs, I)
    ax.axvline(pumpFreq/1.0e9, 0, 1, color = 'r', linestyle = '--')
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    [i.set_linewidth(3.) for i in ax.spines.itervalues()]
    ax.plot(np.round(freqs/1.0e9, 2), gain, linewidth = 3)
    plt.title(r'W-Band Phase shifter Freq Response, f$_{pump}$ = '\
               + str(pumpFreq/1.0e9) + ' GHz', fontsize = 18)
    plt.ylabel('Gain (dB)', fontsize = 18)
    #plt.xlabel(r'$\omega_{s}$ / $\omega_{p}$', fontsize = 18)
    plt.xlabel('frequency (GHz)', fontsize = 18)
    plt.grid()
    plt.savefig('./freqResponse.png', bbox_inches='tight')
    return

############################################
# FTS simultions
###########################################

def interferogram(duration, ramp_freq, samp_freq, sig_freqs, plot = False):
    """Interferogram for one or more freqs, vs DC bias current"""
    n_ramps = np.int(np.ceil(duration * ramp_freq))
    samps_per_ramp = np.int(np.ceil(samp_freq / ramp_freq))
    Ivals = np.linspace(1.0e-10, Istar, samps_per_ramp)
    Irange = []
    for i in range(n_ramps):
        if (i % 2 == 0):
            Irange.append(Ivals)
        else:
            Irange.append(np.flip(Ivals, 0))
    Irange = np.array(Irange)
    Irange = np.hstack(Irange)
    N = len(Irange)
    t = np.linspace(0, duration, N)
    # v = c / np.sqrt(er)
    I0 = 1.
    v = vPhase(Ltot(Irange)/l_max, Cap())
    #delta = l_max*ramp_freq*t
    pshift = np.zeros((len(sig_freqs), N))
    cal_factor = v/(ramp_freq*l_max) # convert spectrum freqs to sig freqs
    I_p = np.zeros((len(sig_freqs), N))
    f_igram = np.zeros((len(sig_freqs), N))
    for i in range(len(sig_freqs)):
        pshift[i] = phi(sig_freqs[i], l_max, v)
        f_igram[i] = sig_freqs[i] * ramp_freq * l_max / v
        #print f_igram
        #pshift[i] = (2*np.pi*sig_freqs[i]*delta) / v
        pshift[i] = pshift[i] - pshift[i][0]
        #print pshift[i]
        I_p[i] = (I0/2.) + (I0/2.)*np.cos(2*np.pi*f_igram[i])
    igram = I_p.sum(axis = 0)
    #f, t, Sxx = spectrogram(igram, samp_freq)
    if plot:
        plt.figure(figsize=(10.24, 7.68), dpi = 100)
        ax = plt.gca()
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        [i.set_linewidth(3.) for i in ax.spines.itervalues()]
        #ax.pcolormesh(t, f, Sxx)
        #for i in range(len(sig_freqs)):
            #ax.plot(t, f_igram[i])
            #ax.plot(t, f_igram[i]*cal_factor)
        #for i in range(len(sig_freqs)):
        #    ax.plot(t, pshift[i])
        #raw_input()
        ax.plot(t, igram)
        #ax.scatter(t, igram)
        ax.set_xlabel('time (s)', fontsize = 18, fontweight='bold')
        ax.set_ylabel('Normalized Intensity', fontsize = 18, fontweight='bold')
        plt.title('Interferogram', fontsize = 18, fontweight = 'bold')
        plt.grid()
    return igram, cal_factor

def spectrum(igram, cal_factor):
    """Plot a spectrum from the interferogram"""
    spec = (2*np.fft.ifft(igram).real[:len(igram)/2])
    freq = np.fft.rfftfreq(igram.size, 1./samp_freq)[1:-1]
    #freq *= cal_factor
    plt.figure(figsize=(10.24, 7.68), dpi = 100)
    ax = plt.gca()
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    [i.set_linewidth(3.) for i in ax.spines.itervalues()]
    ax.plot(freq/1.0e9, spec[1:], linewidth = 2)
    ax.set_xlabel('freq (GHz)', fontsize = 18, fontweight='bold')
    ax.set_ylabel('Power', fontsize = 18, fontweight='bold')
    plt.title('Signal Spectrum', fontsize = 18, fontweight = 'bold')
    plt.grid()
    return

############################
# MAIN
############################
plt.ion()

# plot the phase shift for all freqs, v. bias current
#pShift(freqs, I, plot = True)

# plot the phase shift for all freqs, plus the quadratic gain, v. bias current
#plotPshiftGain(freqs, I)

# plot the parametric gain for a signal at 90 GHz, pump at 90.01 #GHz
#parampGain(Istar, l_max, 90.0e9, 90.01e9)

# plot the paramp's frequency response
#parampFreqResponse(l_max, 90.0e9)

#pShifterFreqResponse(l_max, np.array([0.5*Istar]), 90.0e9)

# plot an interferogram
#samp_freq = 1000 # sampling frequency, Hz
#duration = 1 # s, measurement window
samp_freq = 1000 # sampling frequency, Hz
duration = 2 # s, measurement window
freq_res = (c/(2*l_max)) # Hz
#sig_freqs = np.linspace(75., 112.5, 10)*1.0e9
sig_freqs = np.array([75.0e9])
#ramp_freq = 2. # Hz
ramp_freq = 1.
#I_vals = np.hstack((I, np.flip(I, 0)))
interferogram(duration, ramp_freq, samp_freq, sig_freqs, plot = True)
igram, cal_factor = interferogram(duration, ramp_freq, samp_freq, sig_freqs)

# plot a spectrum from interferogram
spectrum(igram, cal_factor)
