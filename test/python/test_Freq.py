import unittest as ut
import OmegaMaxEnt_TRIQS as OT
from math import ceil, exp, sqrt, pi
from triqs.dos import HilbertTransform, DOSFromFunction
import warnings
from triqs.gf import *
import numpy as np
import os
import shutil as su
warnings.simplefilter(action='ignore', category=FutureWarning)

test_dir_name="test_dir_1"

np.random.seed(1)

tol_int_diffA=0.05

Npts_dos=1000

inter_mode=False
save_figs=False

err=1e-5
err_abs=1e-10
beta=50

R_iw_W=5

W=4
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

def spectr_val(w):
    W = np.sum(wgt)
    v = 0
    for i in range(0,Npks):
        v = v + (wgt[i] / sd[i]) * exp(-(w - cw[i]) * (w - cw[i]) / (2 * sd[i] * sd[i]))

    return v / (W * sqrt(2 * pi))

Nw=int((wr-wl)/dw)+1
w=dw*np.array(list(range(0,Nw)))+wl

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

        errGr=err * np.absolute(G.data.real)
        errGi=err * np.absolute(G.data.imag)

        for i in range(0,2*n_iwn):
            if errGr[i]<err_abs:
                errGr[i] = err_abs

        G.data.real =G.data.real + errGr * np.random.randn(2*n_iwn)
        G.data.imag =G.data.imag + errGi * np.random.randn(2*n_iwn)

        ERRG = errGr + 1j * errGi

        if not os.path.exists(test_dir_name):
            os.mkdir(test_dir_name)
        os.chdir(test_dir_name)

        GR=OT.compute_GfReFreq(G, ERR=ERRG, interactive_mode=inter_mode, save_figures_data=save_figs, output_grid_params=[wl, dw, wr], comp_grid_params=[dw_comp, SW], name="$G_{ME}$")

        os.chdir("..")
        su.rmtree(test_dir_name)

        if isinstance(GR, GfReFreq):
            Aw_me = -GR.data.imag / pi

            int_diffA = dw * sum(np.absolute(Aw_me - Aw))

            print(int_diffA)

            self.assertLess(int_diffA, tol_int_diffA)
        else:
            self.assertTrue(False)

if __name__ == '__main__':
    ut.main()