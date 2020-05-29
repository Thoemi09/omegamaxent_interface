import unittest as ut
import OmegaMaxEnt_TRIQS as OT
from math import ceil, exp, sqrt, pi
from triqs.dos import HilbertTransform, DOSFromFunction
import warnings
from triqs.gf import *
import numpy as np
import os
import shutil as su
# from matplotlib import pyplot as plt
import scipy.integrate as integ
warnings.simplefilter(action='ignore', category=FutureWarning)

test_dir_name="test_dir_8"

delete_test_dir=True
displ_adv_prep_figs=False
preprocess_only=False
inter_mode=False
save_figs=False
compute_Pade=False
eta=1e-5

statistic='Boson'

ref_spectr_file="reference_spectrum.dat"

R_sv_min=1e-30

np.random.seed(0)

tol_int_diffA=0.05

tol_abs=1e-10
tol_rel=1e-6
N_interv_max=1000
# Npts_dos=1000

err=1e-5
beta=50.0

R_iw_W=5

W=8
cw=[-2, 1]
sd=[1, 0.7]
wgt=[1, 1]
Npks=len(cw)

wmin=-10.0
wmax=10.0

wl=-7
wr=7
dw=0.01

dw_comp=0
SW=0
SC=0

wnmax=W*R_iw_W
nmax=int(ceil(beta*wnmax/(2*pi)))
n_iwn=nmax+1

def sum_gaussians(w):
    W = sum(wgt)
    v = 0
    for i in range(0, Npks):
        v = v + (wgt[i] / sd[i]) * exp(-(w - cw[i]) * (w - cw[i]) / (2 * sd[i] * sd[i]))

    return v / (W * sqrt(2 * pi))

def spectr_val(w):
    return w*sum_gaussians(w)

Nw=int((wr-wl)/dw)+1
w=dw*np.array(list(range(0,Nw)))+wl

Aw=np.zeros(Nw)

eps=1e-4

for i in range(0,Nw):
    if abs(w[i])>eps:
        Aw[i]=spectr_val(w[i])/w[i]
    elif w[i]<0:
        Aw[i] = -spectr_val(-eps) / eps
    else:
        Aw[i] = spectr_val(eps) / eps

A_data=np.concatenate((np.array([w]),2*pi*np.array([Aw])),axis=0)
A_data=A_data.transpose()

def integ_sum_gaussians(w,t):
    A = spectr_val(w)
    if w>eps:
        integ=-A*exp(-w*t)/(1-exp(-beta*w))
    elif w<-eps:
        integ=-A*exp(w*(beta-t))/(exp(beta*w)-1)
    else:
        integ=-exp(-w * t)*sum_gaussians(w)/beta
    return integ

def integ_moment(w,n):
    A = spectr_val(w)
    return (w**n)*A

Ntau=n_iwn

dtau=beta/(Ntau-1)

tau=np.arange(0,beta+dtau,dtau)

Gt=np.zeros(Ntau)
erGt=np.zeros(Ntau)

M=np.zeros(4)

for n in range(0,4):
    args=(n,)
    Mtmp=integ.quad(integ_moment, wmin, wmax, args=args, limit=N_interv_max, points=cw, epsabs=tol_abs, epsrel=tol_rel)
    M[n]=Mtmp[0]

print("M="+str(M))

for i in range(0,Ntau):
    args=(tau[i],)
    G_tmp = integ.quad(integ_sum_gaussians, wmin, wmax, args=args, limit=N_interv_max, points=cw, epsabs=tol_abs, epsrel=tol_rel)
    Gt[i]=G_tmp[0]
    erGt[i]=G_tmp[1]

# plt.figure(1)
# plt.plot(tau,Gt)
#
# plt.figure(2)
# plt.plot(tau,erGt)
# plt.show()

class OmegaMaxEnt_test(ut.TestCase):

    def runTest(self):

        # d = DOSFromFunction(spectr_val, wmin, wmax, Npts_dos)
        # Giw = GfImFreq(indices=[0], beta=beta, n_points=n_iwn, statistic=statistic)
        # Sigma0 = GfImFreq(indices=[0], beta=beta, n_points=n_iwn, statistic=statistic)
        # Sigma0.zero()
        # Giw << HilbertTransform(d)(Sigma=Sigma0, mu=0.)
        #
        # Gtau = GfImTime(indices=[0], beta=beta, statistic=statistic)
        # Gtau << InverseFourier(Giw)

        # Gtau = Gtau[0, 0]
        #
        # Gt = Gtau.data.real[::9]

        # Ntau = len(Gt)

        G = GfImTime(target_shape=(), beta=beta, n_points=Ntau, statistic=statistic)
        G.data.real = Gt

        G.data.real = G.data.real + err * np.random.randn(Ntau)

        if not os.path.exists(test_dir_name):
            os.mkdir(test_dir_name)
        os.chdir(test_dir_name)

        np.savetxt(ref_spectr_file, A_data)

        GR=OT.compute_GfReFreq(G,
                               interactive_mode=inter_mode,
                               save_figures_data=save_figs,
                               output_grid_params=[wl, dw, wr],
                               comp_grid_params=[dw_comp, SW],
                               name="$G_{ME}$",
                               displ_adv_preproc_figs=displ_adv_prep_figs,
                               preprocess_only=preprocess_only,
                               ref_spectrum=ref_spectr_file,
                               compute_Pade=compute_Pade,
                               eta_Pade=eta)

        os.chdir("..")
        if delete_test_dir:
            su.rmtree(test_dir_name)

        if isinstance(GR, GfReFreq):
            Aw_me=np.zeros(Nw)

            for i in range(0, Nw):
                if abs(w[i]) > eps:
                    Aw_me[i] = -GR.data[i].imag / (w[i]*pi)
                else:
                    Aw_me[i] = -(GR.data[i+1].imag-GR.data[i-1].imag)/((w[i+1]-w[i-1])*pi)

            int_diffA=dw*sum(np.absolute(Aw_me-Aw))

            print(int_diffA)

            self.assertLess(int_diffA, tol_int_diffA)
        else:
            self.assertTrue(False)

if __name__ == '__main__':
     ut.main()