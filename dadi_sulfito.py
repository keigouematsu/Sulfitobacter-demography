import dadi
from dadi import Numerics, PhiManip, Integration,Misc
from dadi.Spectrum_mod import Spectrum
import random
import pylab
import dadi.cuda

fs = dadi.Spectrum.from_file("./clade12-clade9-clade1-clade2.sfs")
ns = fs.sample_sizes

def split_nomig(params, ns, pts):

    #6 parameters	
    nu1, nuA1, nuA2, nu2, nu3, nu4, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA1, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nuA2, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    # Split off, with contribution entirely from
    # population 3
    phi = PhiManip.phi_3D_to_4D(phi, 0,0, xx,xx,xx,xx)
    
    phi = Integration.four_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, nu4=nu4, m12=0, m13=0, m14=0, m21=0, m23=0, m24=0, m31=0, m32=0, m34=0, m41=0, m42=0, m43=0)


    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx,xx))
    return fs


pts_l = [10]

# Make the extrapolating version of our demographic model function.
dadi.cuda_enabled(True)
func = split_nomig
func_ex = dadi.Numerics.make_extrap_log_func(func)


# initial parameters
p0 = [ 0.01 ,  0.1  ,  0.1  ,  0.01 ,  0.01 ,  0.01 ,  0.002 ,  0.05  ,  0.02]

p_opt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,fixed_params = [None,  None,  None ,  None,  None,  None, None, None, None], verbose=1,epsilon=0.01,  maxiter=10)
print(p_opt)
print('Best-fit parameters: {0}'.format(p_opt))
model = func_ex(p_opt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, fs)
theta = dadi.Inference.optimal_sfs_scaling(model,fs)

print(ll_opt, "\t", p_opt[0], "\t", p_opt[1], "\t", p_opt[2], "\t", p_opt[3], "\t", p_opt[4], "\t", p_opt[5], p_opt[6], "\t", p_opt[7], "\t", p_opt[8], theta)


#2nd round
#-6192.07    , array([ 0.00371792 ,  0.00587973 ,  0.0143727  ,  0.000223558,  0.00441557 ,  0.00838613 ,  0.000335586,  0.00210873 ,  0.0289193 
pts_l = [20]

# Make the extrapolating version of our demographic model function.
dadi.cuda_enabled(True)
func = split_nomig
func_ex = dadi.Numerics.make_extrap_log_func(func)

# initial parameters
p0 = [ 0.00371792 ,  0.00587973 ,  0.0143727  ,  0.000223558,  0.00441557 ,  0.00838613 ,  0.000335586,  0.00210873 ,  0.0289193]
p_opt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,fixed_params = [None ,  None,  None ,  None,None,None  ,  None,None,None], verbose=1,epsilon=0.001,  maxiter=10)
print(p_opt)
print('Best-fit parameters: {0}'.format(p_opt))
model = func_ex(p_opt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, fs)
theta = dadi.Inference.optimal_sfs_scaling(model,fs)

print(ll_opt, "\t", p_opt[0], "\t", p_opt[1], "\t", p_opt[2], "\t", p_opt[3], "\t", p_opt[4], "\t", p_opt[5], p_opt[6], "\t", p_opt[7], "\t", p_opt[8], theta)