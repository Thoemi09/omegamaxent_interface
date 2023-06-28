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

test_dir_name="test_dir_5"

np.random.seed(1)

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

dw=0.01
Npts_dos=int((wmax-wmin)/dw+1)

wl=-15
wr=15
dw=0.005

dw_comp=0
SW=0
SC=0

use_param_grid=True
param_grid=[-4, 0.2, -1.6, 0.06, -0.2, 0.01, 0.2, 0.05, 1.5, 0.1, 3, 0.2, 7]

wnmax=W*R_iw_W

nmax=int(ceil(beta*wnmax/(2*pi)))

ind=np.array(list(range(0,nmax)))
wn=(2*ind+1)*pi/beta

n_iwn=len(wn)

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
        G = GfImFreq(target_shape=[1,1], beta=beta, n_points=n_iwn)
        Sigma0 = GfImFreq(target_shape=[1,1], beta=beta, n_points=n_iwn)
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

        ERRG = errGr + 1j * errGi

        if not os.path.exists(test_dir_name):
            os.mkdir(test_dir_name)
        os.chdir(test_dir_name)

         # GR = OT.compute_GfReFreq(G, ERR=ERRG, interactive_mode=inter_mode, save_figures_data=save_figs, output_grid_params=[wl, dw, wr], comp_grid_params=[dw_comp, SW], non_uniform_grid=nu_grid, name="$G_{ME}$")


        GR = OT.compute_GfReFreq(G, ERR=ERRG, interactive_mode=inter_mode, save_figures_data=save_figs, output_grid_params=[wl, dw, wr],  use_parameterized_grid=use_param_grid, parameterized_grid_params=param_grid, name="$G_{ME}$")

        os.chdir("..")
        su.rmtree(test_dir_name)

        if isinstance(GR, GfReFreq):
            Aw_me=-GR.data.imag/pi

            int_diffA=dw*sum(np.absolute(Aw_me-Aw))

            print(int_diffA)

            self.assertLess(int_diffA, tol_int_diffA)
        else:
            self.assertTrue(False)

if __name__ == '__main__':
    ut.main()
