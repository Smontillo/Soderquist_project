from scipy import integrate, linalg
import json
import numpy as np
import armadillo as arma
from bath_gen_Drude_PSD import generate

# ==============================================================================================
#                                       Global Parameters     
# ==============================================================================================

conv = 27.211397                            # 1 a.u. = 27.211397 eV
fstoau = 41.341                           # 1 fs = 41.341 a.u.
cmtoau = 4.556335e-06                     # 1 cm^-1 = 4.556335e-06 a.u.
autoK = 3.1577464e+05                     # 1 au = 3.1577464e+05 K
kcal_to_au = 1.5936e-03                     # 1 kcal/mol = 1.5936e-3 a.u.
fstocm = 3.33333E4
kB = 3.1668E-6   # au/K
# Model parameters
"""
the model Arkajit used, see orginal paper: JCP 140, 174105 (2014)
"""

def DW(x,m,wDW):
    Eb = 2120 * cmtoau
    V = -(m*wDW**2 / 2) * x**2 + (m**2 * wDW**4 / (16 * Eb)) * x**4
    return V - min(V)

def T(x,m):
    dx = x[1] - x[0]
    N = len(x)
    K = np.pi/dx
    Kin = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            if i == j:
                Kin[i,j] = K**2/3 * (1 + 2/N**2)
            else:
                Kin[i,j] = 2*K**2/N**2 * (-1)**(j-i)/(np.sin(np.pi * (j-i)/N))**2 
    return 1/(2*m) * Kin

def DVR(x,m,wDW):
    V = DW(x,m,wDW)
    V = np.diag(V)
    K = T(x,m)
    E, V = np.linalg.eigh(V+K)
    return E, V

def get_Hs(nDW,lam_Rxn, lam_cav):
    H_DW = np.diag(EDW[:nDW])
    H_bath = (lam_Rxn + lam_cav) * get_Rx2(nDW) 
    return H_DW + H_bath

def get_rho0():
    beta = 1/temp
    rho_DW = np.zeros((nDW,nDW), dtype = 'complex')

    # Reactant well
    rho_DW[0,0] = rho_DW[1,1] = 0.5 + 0j
    rho_DW[1,0] = rho_DW[0,1] = 0.5 + 0j
    
    print(rho_DW)
    return rho_DW

def trace_HO(c):
    rho = c.reshape((nDW*nHO,nDW*nHO))
    rho_t = np.zeros((nDW,nDW), dtype = 'complex')
    for k in range(nDW):
        for l in range(nDW):
            for n in range(nHO):
                rho_t[k,l] += rho[(nDW)*n+k,(nDW)*n+l] 
    return rho_t

def trace_DW(c):
    rho = c.reshape((nDW*nHO,nDW*nHO))
    rho_t = np.zeros((nHO,nHO), dtype = 'complex')
    for k in range(nHO):
        for l in range(nHO):
            for m in range(nDW):
                rho_t[k,l] += rho[(nDW)*k+m,(nDW)*l+m] 
    return rho_t

def pop_DW(c):
    rho_DW = trace_HO(c)
    return np.diag(rho_DW)

def pop_HO(c):
    rho_HO = trace_DW(c)
    return np.diag(rho_HO)

def get_left_overlap(state1, state2):
    P = 0.0
    center = int(len(VDW[:, 0])/2) +1
    dat = VDW[:center, state1].conjugate() * VDW[:center, state2]
    P += np.trapz(dat,x[:center],dx)
    return P

def left_pop(c):
    PRt = 0
    DW = trace_HO(c)
    for i in range(nDW):
        for j in range(nDW):
            PRt += get_left_overlap(i,j) * DW[i,j]
    return PRt

def delta(m, n):
    return 1 if m == n else 0

def get_Rx(nDW):
    pos_DW = np.zeros((nDW,nDW), dtype = 'complex')
    for j in range(nDW):
        for i in range(nDW):
            avg_pos = VDW[:,j].conjugate() * x * VDW[:,i]
            pos_DW[j,i] = np.trapz(avg_pos,x,dx)
    return pos_DW

def get_Rx2(nDW):
    pos_DW = np.zeros((nDW,nDW), dtype = 'complex')
    for j in range(nDW):
        for i in range(nDW):
            avg_pos = VDW[:,j].conjugate() * x**2 * VDW[:,i]
            pos_DW[j,i] = np.trapz(avg_pos,x,dx)
    return pos_DW


# Bath-matter coupling
def get_Qx(nDW):  # R = R_ij |v_i >< v_j|
    return get_Rx(nDW)

# ==============================================================================================
# Create the DW and HO potentials
# ==============================================================================================
nDW = 5
mDW = 1836
wDW = 1000 * cmtoau
N = 1024
L = 3.0
x = np.linspace(-L,L,N)
dx = x[1] - x[0]

EDW,VDW = DVR(x,mDW,wDW)
Normx = np.trapz(VDW[:,0].conjugate() * VDW[:,0],x,dx)
VDW = VDW/(Normx)**0.5
VDW = np.array(VDW, dtype = 'complex')
VDW[:,0] = -VDW[:,0]

temp = 300 / autoK  
# ==============================================================================================
#                                    Summary of parameters     
# ==============================================================================================
class parameters:

    # ===== DEOM propagation scheme =====
    dt = 0.05 * fstoau
    t = 3000 * fstoau   # plateau time as 20ps for HEOM
    nt = int(t / dt)
    nskip = 10

    lmax = 8
    nmax = 1000000
    ferr = 1.0e-07

    # ===== number of system states =====
    x = x
    wDW = wDW
    mDW = mDW
    nDW = nDW
    Nstates = nDW 

    # ===== Drude-Lorentz model =====
    temp = 300 / autoK                 # temperature
    nmod = 2                           # number of dissipation modes (C-Q-R)

    rho0 = get_rho0() #rho_Chen(nDW,nHO) #

    # Bath I (Rxn), Drude-Lorentz model
    gam_Rxn = 200 * cmtoau                      # bath characteristic frequency
    ratio = 0.25                                # the value of etas / omega_b, tune it from 0.02 to 2.0
    lam_Rxn = ratio * mDW * wDW * gam_Rxn/2           # reorganization energy
 
    # Bath II Cav
    wc = 1000 * cmtoau
    eta_c  = 0.05
    tau_c = 1000 * fstoau
    gam_c = 1/tau_c
    lam_c = eta_c**2 * wc

    lam = np.array([lam_Rxn, lam_c])
    gam = np.array([gam_Rxn, gam_c])

    # PSD scheme
    pade    = 1                          # 1 for [N-1/N], 2 for [N/N], 3 for [N+1/N]
    npsd    = 5                      # number of Pade terms

    # ===== Get the subspace information ===== 
    # eigenvalue = eigenvalue
    VDW = VDW 
    V = DW(x,mDW,wDW)

    # ===== Build the bath-free Hamiltonian, dissipation operators, and initial DM in the subspace =====
    # There are three bath ---> DW bath, HO bath and Cavity 'bath'
    # DW --> x, HO --> y

    Qsx = get_Qx(nDW) 
    print(lam)
    print((np.round(np.real(rho0),4)))
    print('Rxn frequency', wDW/cmtoau)
    print('Hamiltonian matrix')
    print(np.round(np.real(get_Hs(nDW,lam_Rxn, lam_c)/cmtoau),1))
    print('Rx')
    print(np.round(np.real(get_Rx(nDW)),2))
    print(EDW[:4]/cmtoau)

# ==============================================================================================
#                                         Main Program     
# ==============================================================================================

if __name__ == '__main__':

    with open('default.json') as f:
        ini = json.load(f)

    # passing parameters
    # bath
    temp = parameters.temp
    nmod = parameters.nmod
    lam = parameters.lam
    gam = parameters.gam
    wc = parameters.wc
    pade = parameters.pade
    npsd = parameters.npsd
    # system
    Nstates = parameters.Nstates
    hams = get_Hs(nDW, lam[0], lam[1])
    rho0 = parameters.rho0
    # system-bath
    Qsx = parameters.Qsx
    # DEOM
    dt = parameters.dt
    nt = parameters.nt
    nskip = parameters.nskip
    lmax = parameters.lmax
    nmax = parameters.nmax
    ferr = parameters.ferr

# ==============================================================================================================================
    # hidx
    ini['hidx']['trun'] = 0
    ini['hidx']['lmax'] = lmax
    ini['hidx']['nmax'] = nmax
    ini['hidx']['ferr'] = ferr

	# bath PSD
    ini['bath']['temp'] = temp
    ini['bath']['nmod'] = nmod
    ini['bath']['pade'] = pade
    ini['bath']['npsd'] = npsd
    ini['bath']['jomg'] = [{"jdru":[(lam[i], gam[i])]} for i in range(nmod-1)] 
                                                                               
    jomg = ini['bath']['jomg']
    nind = 0
    for m in range(nmod):       # one mode is treated by PFD
        try:
            ndru = len(jomg[m]['jdru'])
        except:
            ndru = 0
        try:
            nsdr = len(jomg[m]['jsdr'])
        except:
            nsdr = 0
        nper = ndru + 2 * nsdr + npsd
        nind += nper
                                                                               
    etal_1, etar_1, etaa_1, expn_1, delr_1 = generate (temp, npsd, pade, jomg)

    # bath II with PSD
    ini['bath']['temp'] = temp                                                  
    ini['bath']['nmod'] = nmod
    ini['bath']['pade'] = pade
    ini['bath']['npsd'] = npsd + 1
    ini['bath']['jomg'] = [{"jsdr":[(lam[1], wc, gam[1])]}] 

    jomg = ini['bath']['jomg']

    etal_2, etar_2, etaa_2, expn_2, delr_2 = generate (temp, npsd + 1, pade, jomg)

    mode = np.zeros((nind + 2 + npsd + 1), dtype = int)
    for i in range(nind, nind + npsd + 3):
        mode[i] = 1

    delr = np.append(delr_1, delr_2)
    etal = np.append(etal_1, etal_2)
    etar = np.append(etar_1, etar_2)
    etaa = np.append(etaa_1, etaa_2)
    expn = np.append(expn_1, expn_2)

    arma.arma_write(mode, 'inp_mode.mat')
    arma.arma_write(delr, 'inp_delr.mat')
    arma.arma_write(etal, 'inp_etal.mat')
    arma.arma_write(etar, 'inp_etar.mat')
    arma.arma_write(etaa, 'inp_etaa.mat')
    arma.arma_write(expn, 'inp_expn.mat')

    # two dissipation modes
    qmds = np.zeros((nmod, Nstates, Nstates), dtype = complex)
    qmds[0,:,:] = Qsx * 1                    # the electron-phonon interaction (DW)
    qmds[1,:,:] = Qsx * 1

    arma.arma_write (hams,ini['syst']['hamsFile'])
    arma.arma_write (qmds,ini['syst']['qmdsFile'])
    arma.arma_write (rho0,'inp_rho0.mat')

    jsonInit = {"deom":ini,
                "rhot":{
                    "dt": dt,
                    "nt": nt,
                    "nk": nskip,
					"xpflag": 1,
					"staticErr": 0,
                    "rho0File": "inp_rho0.mat",
                    "sdipFile": "inp_sdip.mat",
                    "pdipFile": "inp_pdip.mat",
					"bdipFile": "inp_bdip.mat"
                },
            }

# ==============================================================================================================================
# ==============================================================================================================================

    # dipoles
    sdip = np.zeros((2,2),dtype=float)
    arma.arma_write(sdip,'inp_sdip.mat')

    pdip = np.zeros((nmod,2,2),dtype=float)
    pdip[0,0,1] = pdip[0,1,0] = 1.0
    arma.arma_write(pdip,'inp_pdip.mat')

    bdip = np.zeros(3,dtype=complex)
#    bdip[0]=-complex(5.00000000e-01,8.66025404e-01)
#    bdip[1]=-complex(5.00000000e-01,-8.66025404e-01)
#    bdip[2]=-complex(7.74596669e+00,0.00000000e+00)
    arma.arma_write(bdip,'inp_bdip.mat')

    with open('input.json','w') as f:
        json.dump(jsonInit,f,indent=4) 
