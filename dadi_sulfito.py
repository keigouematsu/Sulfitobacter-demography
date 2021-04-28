import dadi
from dadi import Numerics, PhiManip, Integration,Misc
from dadi.Spectrum_mod import Spectrum
import random

fs = dadi.Spectrum.from_file("./clade1-clade2-clade3.sfs")
ns = fs.sample_sizes

def split_nomig(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #6 parameters	
    nu1, nuA, nu2, nu3, T1, T2 = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


# the number of grid points for extrapolation.
pts_l = [20]

# Make the extrapolating version of our demographic model function.
func = split_nomig
func_ex = dadi.Numerics.make_extrap_log_func(func)


for i in range (3):
	print('loop round', repr(i))
	# initial parameters
	p0 = [2.18949e-05,  0.000149668,  0.00431911 ,  1.91889e-05,  8.47597e-06,  0.0316052]

	p_opt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, verbose=1, maxiter=10, epsilon=0.01)
	print(p_opt)
	print('Best-fit parameters: {0}'.format(p_opt))
	model = func_ex(p_opt, ns, pts_l)
	ll_opt = dadi.Inference.ll_multinom(model, fs)
	theta = dadi.Inference.optimal_sfs_scaling(model,fs)

	print(ll_opt, "\t", p_opt[0], "\t", p_opt[1], "\t", p_opt[2], "\t", p_opt[3], "\t", p_opt[4], "\t", p_opt[5], theta)
