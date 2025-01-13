from scipy import integrate, linalg
import json
import numpy as np
import armadillo as arma
from bath_gen_Drude_PSD import generate

# ==============================================================================================
#  PHYSICAL CONSTANTS 
# ==============================================================================================
autoeV = 27.211397         # 1 a.u. = 27.211397 eV
fstoau = 41.341            # 1 fs = 41.341 a.u.
cmtoau = 4.556335e-06      # 1 cm^-1 = 4.556335e-06 a.u.
autoK  = 3.1577464e+05     # 1 au = 3.1577464e+05 K
fstocm = 3.33333E4
kB = 3.1668E-6   # au/K

# ==============================================================================================
#  FUNCTIONS
# ==============================================================================================


# BATH - MATTER COUPLING
def get_Qx():  
    return 

# ==============================================================================================
# CREATE THE DOUBLE-WELL EIGENFUNCTIONS AND ENERGIES (AS A FUNCTION OF THE NUMBER OF STATES)
# ==============================================================================================



# ==============================================================================================
#                                    Summary of parameters     
# ==============================================================================================
class parameters:
    # ===== DEOM propagation scheme =====
    dt      = 0.1 * fstoau                     # time step
    t       = 5000 * fstoau                    # final time
    nt      = int(t / dt)                      # number of steps
    nskip   = 10                               # number of steps between saved data

    lmax    = 8                                # maximum HEOM depth
    nmax    = 1000000                          # maximum HEOM functions
    ferr    = 1.0e-07                          # integration error
    # ===== system parameters =====
    x       =                                  # Position operator in position basis
    wDW     =                                  # Top of the barrier frequency
    mDW     =                                  # Mass of the particle
    Nstates =                                  # Number of vibrational states
    # ===== Drude-Lorentz model =====
    temp    = 300 / autoK                      # temperature
    nmod    = 2                                # number of dissipation modes 
    # Bath I (Drude-Lorentz model)
    gam_Rxn = 200 * cmtoau                     # bath characteristic frequency
    ratio   = 0.1                              # the value of etas / omega_b, tune it from 0.02 to 2.0
    lam_Rxn = ratio * mDW * wDW * gam_Rxn/2    # reorganization energy
    # Bath II Cav (Brownian Spectral Density model)
    wc      = 1000 * cmtoau                    # cavity frequency
    eta_c   = 0.05                             # light - matter coupling
    tau_c   = 1000 * fstoau                    # cavity lifetime
    gam_c   = 1/tau_c                          # bath characteristic frequency
    lam_c   = eta_c**2 * wc                    # bath reorganization energy

    lam     = np.array([lam_Rxn, lam_c])       # reorganization energies
    gam     = np.array([gam_Rxn, gam_c])       # characteristic frequencies

    # PSD scheme (Decomposition scheme for the spectral density)
    pade    = 1                                # 1 for [N-1/N], 2 for [N/N], 3 for [N+1/N]
    npsd    = 3                                # number of Pade terms
    # ===== Build the bath-free Hamiltonian, dissipation operators, and initial DM in the subspace =====
    # There are two baths ---> DW bath and Cavity 'bath'
    # DW --> x, HO --> y
    rho0    =                                  # Initial Density Matrix
    Qsx     =                                  # System bath coupling operator
    ham     =                                  # Bath-free Hamiltonian

# ==============================================================================================
#                                         Main Program     
# ==============================================================================================

if __name__ == '__main__':

    with open('default.json') as f:
        ini = json.load(f)

    # passing parameters
    # bath
    temp    = parameters.temp
    nmod    = parameters.nmod
    lam     = parameters.lam
    gam     = parameters.gam
    wc      = parameters.wc
    pade    = parameters.pade
    npsd    = parameters.npsd
    # system
    Nstates = parameters.Nstates
    hams    = parameters.ham
    rho0    = parameters.rho0
    # system-bath
    Qsx     = parameters.Qsx
    # DEOM
    dt      = parameters.dt
    nt      = parameters.nt
    nskip   = parameters.nskip
    lmax    = parameters.lmax
    nmax    = parameters.nmax
    ferr    = parameters.ferr

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
    qmds[0,:,:] = Qsx * 1  # the electron-phonon interaction (DW)
    qmds[1,:,:] = Qsx * 1  # light - matter interaction

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
