import unittest as ut
import OmegaMaxEnt_TRIQS as OT
from math import cos, sin, ceil, exp, sqrt, pi
from pytriqs.dos import HilbertTransform, DOSFromFunction
import warnings
from pytriqs.gf import *
import numpy as np
import os
import shutil as su
warnings.simplefilter(action='ignore', category=FutureWarning)

test_dir_name="test_dir_4"

np.random.seed(0)

tol_int_diffA=0.05

Npts_dos=1000

theta = pi / 8

inter_mode=False
save_figs=False
inv_sym=True

err=1e-5
beta=50

R_iw_W=5

W=4
cw1=[-2, 1]
sd1=[1, 0.7]
wgt1=[1, 1]
Npks1=len(cw1)

cw2=[-1.5, 0, 2]
sd2=[0.8, 0.5, 1.2]
wgt2=[1, 1, 1]
Npks2=len(cw2)

cw3=[-1, 0.2, 1.7]
sd3=[0.6, 0.5, 1]
wgt3=[1, 1, 1]
Npks3=len(cw3)

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

ind=np.array(range(0,nmax+1))
wn=(2*ind+1)*pi/beta

n_iwn=nmax+1

def spectr_val1(w):
    W = np.sum(wgt1)
    v = 0
    for i in range(0,Npks1):
        v = v + (wgt1[i] / sd1[i]) * exp(-(w - cw1[i]) * (w - cw1[i]) / (2 * sd1[i] * sd1[i]))

    return v / (W * sqrt(2 * pi))

def spectr_val2(w):
    W = np.sum(wgt2)
    v = 0
    for i in range(0,Npks2):
        v = v + (wgt2[i] / sd2[i]) * exp(-(w - cw2[i]) * (w - cw2[i]) / (2 * sd2[i] * sd2[i]))

    return v / (W * sqrt(2 * pi))

def spectr_val3(w):
    W = np.sum(wgt3)
    v = 0
    for i in range(0,Npks3):
        v = v + (wgt3[i] / sd3[i]) * exp(-(w - cw3[i]) * (w - cw3[i]) / (2 * sd3[i] * sd3[i]))

    return v / (W * sqrt(2 * pi))

Nw=int((wr-wl)/dw)+1
w=dw*np.array(range(0,Nw))+wl

Aw1=np.zeros(Nw)
Aw2=np.zeros(Nw)
Aw3=np.zeros(Nw)

for i in range(0,Nw):
    Aw1[i] = spectr_val1(w[i])
    Aw2[i] = spectr_val2(w[i])
    Aw3[i] = spectr_val3(w[i])

R = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
Rt=R.transpose()

A00=R[0,0]*Aw1*Rt[0,0]+R[0,1]*Aw2*Rt[1,0]
A01=R[0,0]*Aw1*Rt[0,1]+R[0,1]*Aw2*Rt[1,1]
A11=R[1,0]*Aw1*Rt[0,1]+R[1,1]*Aw2*Rt[1,1]

class OmegaMaxEnt_test_with_error(ut.TestCase):

    def runTest(self):
        Gbl1 = GfImFreq(indices=[0, 1], beta=beta, n_points=n_iwn)

        d1 = DOSFromFunction(spectr_val1, wmin, wmax, Npts_dos)
        G1 = GfImFreq(indices=[0], beta=beta, n_points=n_iwn)
        Sigma0 = GfImFreq(indices=[0], beta=beta, n_points=n_iwn)
        Sigma0.zero()
        G1 << HilbertTransform(d1)(Sigma=Sigma0, mu=0.)

        d2 = DOSFromFunction(spectr_val2, wmin, wmax, Npts_dos)
        G2 = GfImFreq(indices=[0], beta=beta, n_points=n_iwn)
        G2 << HilbertTransform(d2)(Sigma=Sigma0, mu=0.)

        Gbl1[0, 0] = G1[0, 0]
        Gbl1[1, 1] = G2[0, 0]

        Gbl2 = GfImFreq(indices=[0], beta=beta, n_points=n_iwn)

        d3 = DOSFromFunction(spectr_val3, wmin, wmax, Npts_dos)
        Gbl2 << HilbertTransform(d3)(Sigma=Sigma0, mu=0.)

        G_rot = GfImFreq(indices=[0, 1], beta=beta, n_points=n_iwn)
        G_rot.from_L_G_R(R, Gbl1, Rt)

        G_rot.data.real = G_rot.data.real + err * np.reshape(np.random.randn(np.size(G_rot.data.real)),np.shape(G_rot.data.real))
        G_rot.data.imag = G_rot.data.imag + err * np.reshape(np.random.randn(np.size(G_rot.data.real)),np.shape(G_rot.data.real))

        Gbl2.data.real = Gbl2.data.real + err * np.reshape(np.random.randn(np.size(Gbl2.data.real)),np.shape(Gbl2.data.real))
        Gbl2.data.imag = Gbl2.data.imag + err * np.reshape(np.random.randn(np.size(Gbl2.data.real)),np.shape(Gbl2.data.real))

        ERRG = 0

        G = BlockGf(name_list=['bl1', 'bl2'], block_list=[G_rot, Gbl2])

        if not os.path.exists(test_dir_name):
            os.mkdir(test_dir_name)
        os.chdir(test_dir_name)

        GR = OT.compute_GfReFreq(G, ERR=ERRG, interactive_mode=inter_mode, save_figures_data=save_figs, inv_sym=inv_sym,
                                 output_grid_params=[wl, dw, wr], comp_grid_params=[dw_comp, SW], name="$G_{ME}$")

        os.chdir("..")
        su.rmtree(test_dir_name)

        A00_me = -GR['bl1'][0, 0].data.imag / pi
        A01_me = -GR['bl1'][0, 1].data.imag / pi
        A11_me = -GR['bl1'][1, 1].data.imag / pi
        A3_me = -GR['bl2'][0, 0].data.imag / pi

        int_diff_A00=dw*sum(np.absolute(A00_me-A00))
        int_diff_A01=dw*sum(np.absolute(A01_me-A01))
        int_diff_A11=dw*sum(np.absolute(A11_me-A11))
        int_diff_A3=dw * sum(np.absolute(A3_me - Aw3))

        print int_diff_A00
        print int_diff_A01
        print int_diff_A11
        print int_diff_A3

        t00 = int_diff_A00 < tol_int_diffA
        t01 = int_diff_A01 < tol_int_diffA
        t11 = int_diff_A11 < tol_int_diffA
        t3 = int_diff_A3 < tol_int_diffA

        self.assertTrue(t00 and t01 and t11 and t3)

if __name__ == '__main__':
    ut.main()