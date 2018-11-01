import unittest as ut
import OmegaMaxEnt_TRIQS as OT
# import scipy.integrate as integ
from math import ceil, exp, sqrt, pi
from pytriqs.dos import HilbertTransform, DOSFromFunction
# from pytriqs.plot.mpl_interface import oplot
from matplotlib import pyplot as plt
import warnings
from pytriqs.gf import *
import numpy as np
import os
import shutil as su
warnings.simplefilter(action='ignore', category=FutureWarning)

np.random.seed(0)

tol_int_diffA=0.05

inter_mode=False
save_figs=False
nu_grid=True

tol_abs=1e-10
tol_rel=1e-8

err=1e-6
err_abs=1e-10
beta=100

R_iw_W=5

W=8
cw=[-3, -1, 0, 4]
sd=[1, 0.3 ,0.05, 2.5]
wgt=[1, 1, 1, 1]

Npks=len(cw)

wmin=-20.0
wmax=20.0

N_interv_max=10000
Npts_dos=20000

wl=-15
wr=15
dw=0.005

dw_comp=0
SW=0
SC=0

wnmax=W*R_iw_W

nmax=int(ceil(beta*wnmax/(2*pi)))

ind=np.array(range(0,nmax))
wn=(2*ind+1)*pi/beta

n_iwn=len(wn)

def spectr_val(w):
    W = np.sum(wgt)
    v = 0
    for i in range(0,Npks):
        v = v + (wgt[i] / sd[i]) * exp(-(w - cw[i]) * (w - cw[i]) / (2 * sd[i] * sd[i]))

    return v / (W * sqrt(2 * pi))

Nw=int((wr-wl)/dw)+1
w=dw*np.array(range(0,Nw))+wl

Aw=np.zeros(Nw)

for i in range(0,Nw):
    Aw[i]=spectr_val(w[i])

class OmegaMaxEnt_test_with_error(ut.TestCase):

    def runTest(self):

        d = DOSFromFunction(spectr_val, wmin, wmax, Npts_dos)
        G = GfImFreq(indices=[0], beta=beta, n_points=n_iwn)
        Sigma0 = GfImFreq(indices=[0], beta=beta, n_points=n_iwn)
        Sigma0.zero()
        G << HilbertTransform(d)(Sigma=Sigma0, mu=0.)

        G = G[0, 0]

        Gr = G.data.real
        Gi = G.data.imag

        G = GfImFreq(target_shape=(), beta=beta, n_points=n_iwn)

        errGr = err * np.absolute(Gr)
        errGi = err * np.absolute(Gi)

        for i in range(0, 2 * n_iwn):
            if errGr[i] < err_abs:
                errGr[i] = err_abs

        G.data.real = Gr + errGr * np.random.randn(2 * n_iwn)
        G.data.imag = Gi + errGi * np.random.randn(2 * n_iwn)

        errGr = np.array([errGr])
        errGi = np.array([errGi])

        errGr = errGr.transpose()
        errGi = errGi.transpose()

        ERRG = np.concatenate((errGr, errGi), axis=1)

        if not os.path.exists("test_dir"):
            os.mkdir("test_dir")
        os.chdir("test_dir")

        GR = OT.compute_GfReFreq(G, ERR=ERRG, interactive_mode=inter_mode, save_figures_data=save_figs,
                                 grid_params=[wl, dw, wr], comp_grid_params=[dw_comp, SW], non_uniform_grid=nu_grid, name="$G_{ME}$")

        os.chdir("..")
        su.rmtree("test_dir")

        Aw_me=-GR.data.imag/pi

        int_diffA=dw*sum(np.absolute(Aw_me-Aw))

        print int_diffA

        self.assertLess(int_diffA, tol_int_diffA)

if __name__ == '__main__':
    ut.main()