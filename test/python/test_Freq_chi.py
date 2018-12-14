import unittest as ut
import OmegaMaxEnt_TRIQS as OT
from math import ceil, exp, sqrt, pi
from pytriqs.dos import HilbertTransform, DOSFromFunction
import warnings
from pytriqs.gf import *
import numpy as np
import os
import shutil as su
from matplotlib import pyplot as plt
import scipy.integrate as integ
from pytriqs.archive import HDFArchive as HA
warnings.simplefilter(action='ignore', category=FutureWarning)

test_dir_name="test_dir_7"

ref_spectr_file="reference_spectrum.dat"

delete_test_dir=True
inv_sym_time=True
preprocess_only=False
inter_mode=False
save_figs=False
displ_preproc_figs=False
compute_Pade=False
eta=1e-5

statistic='Boson'

np.random.seed(0)

tol_int_diffA=0.05

eps=1e-4

Npts_dos=2000

tol_abs=1e-10
tol_rel=1e-7
N_interv_max=2000

err=1e-5
err_abs=1e-10
beta=50.0

R_iw_W=5

W=8
cw=[-2, -1, 1, 2]
sd=[1, 0.7, 0.7, 1]
wgt=[1, 1, 1, 1]
Npks=len(cw)

wmin=-10.0
wmax=10.0

wl=-7
wr=7
dw=0.01

dw_comp=0
SW=0

wnmax=W*R_iw_W
nmax=int(ceil(beta*wnmax/(2*pi)))
n_iwn=nmax+1

n=np.array(range(0,n_iwn))

wn=2*n*pi/beta

def sum_gaussians(w):
    W = sum(wgt)
    v = 0
    for i in range(0, Npks):
        v = v + (wgt[i] / sd[i]) * exp(-(w - cw[i]) * (w - cw[i]) / (2 * sd[i] * sd[i]))

    return v / (W * sqrt(2 * pi))

def spectr_val(w):
    return w*sum_gaussians(w)

def integ_spectr_re(w,wn):
    if wn:
        integ=-w*w*sum_gaussians(w)/(wn*wn+w*w)
    else:
        integ=-sum_gaussians(w)
    return integ

def integ_spectr_im(w,wn):
    if wn:
        integ=-wn*w*sum_gaussians(w)/(wn*wn+w*w)
    else:
        integ=0
    return integ

Nw=int((wr-wl)/dw)+1
w=dw*np.array(range(0,Nw))+wl

Aw=np.zeros(Nw)

for i in range(0,Nw):
    Aw[i]=sum_gaussians(w[i])

A_data=np.concatenate((np.array([w]),pi*np.array([Aw])),axis=0)
A_data=A_data.transpose()

def integ_moment(w,p):
    A = spectr_val(w)
    return (w**p)*A

M=np.zeros(4)

for p in range(0,4):
    args=(p,)
    Mtmp=integ.quad(integ_moment, wmin, wmax, args=args, limit=N_interv_max, points=cw, epsabs=tol_abs, epsrel=tol_rel)
    M[p]=Mtmp[0]

print "M="
print M


class OmegaMaxEnt_test(ut.TestCase):

    def runTest(self):

        Gr=np.zeros(n_iwn)
        Gi=np.zeros(n_iwn)

        for i in n:
            args=(wn[i],)
            Gtmp=integ.quad(integ_spectr_re, wmin, wmax, args=args, limit=N_interv_max, points=cw, epsabs=tol_abs,
                            epsrel=tol_rel)
            Gr[i]=Gtmp[0]
            Gtmp = integ.quad(integ_spectr_im, wmin, wmax, args=args, limit=N_interv_max, points=cw, epsabs=tol_abs,
                              epsrel=tol_rel)
            Gi[i] = Gtmp[0]

        Gr_n=np.flipud(Gr[1:])
        Gi_n=np.flipud(-Gi[1:])
        Gr=np.concatenate((Gr_n,Gr))
        Gi=np.concatenate((Gi_n,Gi))

        G = GfImFreq(target_shape=(), beta=beta, n_points=n_iwn, statistic=statistic)
        G.data.real=Gr
        G.data.imag=Gi

        # d = DOSFromFunction(spectr_val, wmin, wmax, Npts_dos)
        # G = GfImFreq(indices=[0], beta=beta, n_points=n_iwn, statistic=statistic)
        # Sigma0 = GfImFreq(indices=[0], beta=beta, n_points=n_iwn, statistic=statistic)
        # Sigma0.zero()
        # G << HilbertTransform(d)(Sigma=Sigma0, mu=0.)
        #
        # G=-G[0,0]/2

        # plt.figure(1)
        # plt.plot(G.data.real,'bo')
        # plt.plot(G.data.imag, 'ro')
        # plt.show()

        n_iwn_all=len(G.mesh)

        errGr = err * np.absolute(G.data.real)
        errGi = err * np.absolute(G.data.imag)

        for i in range(0, n_iwn_all):
            if errGi[i] < err_abs:
                errGi[i] = err_abs

        G.data.real = G.data.real + errGr * np.random.randn(n_iwn_all)
        G.data.imag = G.data.imag + errGi * np.random.randn(n_iwn_all)

        ERRG = errGr + 1j * errGi

        if not os.path.exists(test_dir_name):
            os.mkdir(test_dir_name)
        os.chdir(test_dir_name)

        np.savetxt(ref_spectr_file, A_data)

        GR=OT.compute_GfReFreq(G,
                               ERR=ERRG,
                               interactive_mode=inter_mode,
                               save_figures_data=save_figs,
                               output_grid_params=[wl, dw, wr],
                               comp_grid_params=[dw_comp, SW],
                               name="$G_{ME}$",
                               inv_sym_time=inv_sym_time,
                               displ_preproc_figs=displ_preproc_figs,
                               compute_Pade = compute_Pade,
                               eta_Pade = eta,
                               preprocess_only=preprocess_only,
                               ref_spectrum=ref_spectr_file)

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

            print int_diffA

            self.assertLess(int_diffA, tol_int_diffA)
        else:
            self.assertTrue(False)

if __name__ == '__main__':
     ut.main()