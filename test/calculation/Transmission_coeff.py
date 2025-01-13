import numpy as np
import matplotlib.pyplot as plt
from gen_input import parameters, DW

# ==============================================================================================
#                                       Global Parameters     
# ==============================================================================================

conv = 27.211397                            # 1 a.u. = 27.211397 eV
fstoau = 41.341                           # 1 fs = 41.341 a.u.
cmtoau = 4.556335e-06                     # 1 cm^-1 = 4.556335e-06 a.u.
autoK = 3.1577464e+05                     # 1 au = 3.1577464e+05 K
kcaltoau = 1.5936e-03                     # 1 kcal/mol = 1.5936e-3 a.u.
# ==============================================================================================

x = parameters.x
m = parameters.mDW
wDW = parameters.wDW
NStates = parameters.Nstates
VDW = parameters.VDW
nDW = parameters.nDW
V = parameters.V


def rho_t(rho):
    rho_t = np.zeros((points,nDW,nDW), dtype = 'complex')
    for k in range(points):
        for i in range(NStates):
            for j in range(NStates):
                rho_t[k,i,j] = rho[k,NStates * i+j+1]
    return rho_t

def trace_HO(rho):
    rho_t = np.zeros((points,nDW,nDW), dtype = 'complex')
    for j in range(points):
        for k in range(nDW):
            for l in range(nDW):
                for n in range(nHO):
                    rho_t[j,k,l] += rho[j,((nDW)*n+k) * Nstates + (nDW)*n+l + 1] 
    return rho_t

def get_left_overlap(state1, state2):
    P = 0.0
    dx = x[1] - x[0]
    center = int(len(VDW[:, 0])/2) +1
    dat = VDW[:center, state1].conjugate() * VDW[:center, state2]
    P += np.trapz(dat,x[:center],dx)
    return P

def get_right_overlap(state1, state2):
    P = 0.0
    dx = x[1] - x[0]
    center = int(len(VDW[:, 0])/2) +1
    dat = VDW[center:, state1].conjugate() * VDW[center:, state2]
    P += np.trapz(dat,x[center:],dx)
    return P

def left_pop(DW):
    PRt = np.zeros(points, dtype = 'complex')
    for i in range(nDW):
        for j in range(nDW):
            PRt[:] += get_left_overlap(i,j) * DW[:,i,j]
    return PRt

def pop(DW):
    PRt = np.zeros((points,nDW), dtype = 'complex')
    num = range(4)
    for i in num[:2]:
        for j in num[:2]:
            # PRt[:,0] += get_left_overlap(i,j) * DW[:,i,j]
            PRt[:,1] += get_right_overlap(i,j) * DW[:,i,j]

    for i in num[2:]:
        for j in num[2:]:
            PRt[:,2] += get_left_overlap(i,j) * DW[:,i,j]
            PRt[:,3] += get_right_overlap(i,j) * DW[:,i,j]
    return PRt

def rate(DW, time):
    dt = time[1] - time[0]
    PRt = left_pop(DW)
    rate = np.zeros((len(PRt) - 1), dtype = float)
    for n in range(1, len(rate)):
        dPr = ((PRt[n] - PRt[n-1])/dt)
        rate[n] = np.real(dPr/(1.0 - 2.0*PRt[n]))
    return rate


plt.plot(x,VDW[:,0:5])
plt.savefig('images/wf.png')
plt.close()  

freq = ['0.1']#['0.05', '0.1', '0.25', '0.5', '1.0', '1.5', '2.0', '2.5']
kappas = np.zeros(len(freq))
Plot = True

for n in range(len(freq)):
    rhot = np.loadtxt(f'summary/{freq[n]}.dat')
    time = rhot[:,0]
    points = len(time)
    DWt = rho_t(rhot)
    
    rates = rate(DWt, time)
    kappas[n] = rates[-1]
    if Plot:
        plt.plot(rhot[1:,0]/fstoau,rates*fstoau, label = f'{freq[n]}')
        plt.ylabel('$\kappa$')
        plt.xlabel('Steps')
plt.legend()
plt.savefig(f'images/kappa.png', dpi=300, bbox_inches='tight')
plt.close()

# comp = np.loadtxt('compare.txt')
# plt.plot(comp[:,0],comp[:,1], c = 'blue')#, label = 'Lindoy (Digitized)')
w = np.zeros(len(freq))
for i in range(len(freq)):
    w[i] = float(freq[i])
plt.scatter(w,kappas, c = 'red', label = 'Implementation')

plt.ylabel('$\kappa$')
plt.xlabel('$ω_Q$')
# plt.xlim(800,1500)
plt.legend()
plt.savefig('Rates.png', dpi=300, bbox_inches='tight')
plt.close()

for n in range(len(freq)):
    rhot = np.loadtxt(f'summary/{freq[n]}.dat')
    time = rhot[:,0]
    points = len(rhot[:,0])
    
    DWt = rho_t(rhot)
    L_pop = left_pop(DWt)

    if Plot:
        plt.plot(rhot[:,0]/fstoau,L_pop, label = f'{freq[n]}')
        plt.ylabel('$\kappa$')
        plt.xlabel('Time (fs)')
plt.legend()
plt.savefig(f'images/Left.png', dpi=300, bbox_inches='tight')
plt.close()

for n in range(len(freq)):
    rhot = np.loadtxt(f'summary/{freq[n]}.dat')
    time = rhot[:,0]
    points = len(rhot[:,0])
    
    DWt = rho_t(rhot)
    L_pop = left_pop(DWt)
    dia = pop(DWt)

    for e in range(5):
        plt.plot(time/fstoau,dia[:,e+1], label = f'State {e+1}')
    plt.xlabel('Time (fs)')
    plt.ylabel('Population')
    plt.title('Double well diabats population')
    plt.legend()
    plt.savefig(f'images/dia_{freq[n]}.png', dpi=300, bbox_inches='tight')
    plt.close()
    pop_DWt = np.zeros((points,nDW+1))
    for h in range(len(time)):
        tDW = (np.diag(DWt[h,:,:]))
        pop_DWt[h,:-1] = tDW
        pop_DWt[h,-1] = sum(tDW)

    plt.plot(time/fstoau,pop_DWt)
    plt.savefig(f'images/DW_pop_{freq[n]}.png', dpi=300, bbox_inches='tight')
    plt.close()


