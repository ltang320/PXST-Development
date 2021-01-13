import os,sys
import numpy as np
import yaml
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math

"""
------------
import Parameters

"""
with open('configure.yml','r') as conf_para:
    conf_para = yaml.load(conf_para,Loader=yaml.FullLoader)

"""
wavefront_initialize

input
------------
pixelsize_x = 55e-06,pixelsize_y=55e-06,
fs_size = 2000,ss_size = 2000,
focus_x = 1.2e-3,focus_y = 1.0e-3,
defocus = 400e-6,
det_dist = 14e-03,
ap_x = 40e-06, ap_y= 40e-6,
wl = 7.29e-11,
amplitude_value=0.0

output
------------
x_arr,y_arr,xx_arr,yy_arr,wf_dec

"""
def wavefront_initialize(pixelsize_x = 55e-06,pixelsize_y = 55e-06,fs_size = 2000,ss_size = 2000,focus_x = 1.2e-3,focus_y = 1.0e-3,defocus = 400e-6, det_dist = 14e-03, ap_x = 40e-06, ap_y= 40e-6,wl = 7.29e-11,amplitude_value=0.0):

    wf_dec = np.zeros((ss_size,fs_size),dtype='complex')
    wf_dec += amplitude_value


    # the range of detector plane(x-axis,y-axis)
    xx_span = fs_size * pixelsize_x
    yy_span = ss_size * pixelsize_y

    # the range of object plane(x-axis,y-axis)
    x_span = 1.6 * ap_x / focus_x * defocus
    y_span = 1.6 * ap_y / focus_y * defocus
    # the sample rate in the object plane
    n_x = int(x_span * xx_span / wl / det_dist)
    n_y = int(y_span * yy_span / wl / det_dist)

    # Initializing coordinate arrays
    # coordinate in object plane
    x_arr = np.linspace(-x_span / 2, x_span / 2, n_x)
    y_arr = np.linspace(-y_span / 2, y_span / 2, n_y)

    wf_ini = np.zeros((n_x,n_y),dtype='complex')

    # coordinate in detector plan
    xx_arr = np.linspace(-xx_span / 2, xx_span / 2, fs_size, endpoint=False)
    yy_arr = np.linspace(-yy_span / 2, yy_span / 2, ss_size, endpoint=False)

    return x_arr,y_arr,xx_arr,yy_arr,wf_ini,wf_dec

def aperture_type(x_arr,y_arr,wf_ini,ap_x = 40e-06,ap_y = 40e-06,type =0,sigma = 2.35,focus_x = 1.2e-3,focus_y = 1.2e-3, defocus = 400e-06):
    ap_xx = x_arr[:,np.newaxis]
    ap_yy = y_arr[np.newaxis,:]
    diameter_x = ap_x/focus_x*defocus
    diameter_y = ap_y/focus_y*defocus
    filter = np.zeros_like(wf_ini)
    if type == 0: # Rectangular aperture_type
        radius_x = (diameter_x/2)
        radius_y = (diameter_y/2)
        filter_illuminated_indices = np.where((np.abs(ap_xx)<radius_x) & (np.abs(ap_yy)<radius_y))
        filter[filter_illuminated_indices] = 1.0
    elif type == 1: # Circular aperture_type
        radius = np.sqrt((diameter_x/2)**2 + (diameter_y/2)**2)
        filter_illuminated_indices = np.where((ap_xx**2 + ap_yy**2)<radius**2)
        filter[filter_illuminated_indices] = 1.0
    elif type == 2: # Gaussian_type
        sigma_x = diameter_x/2.35
        sigma_y = diameter_y/2.35
        filter = np.sqrt(np.exp(-pow(ap_xx/sigma_x,2)/2-pow((ap_yy/sigma_y),2)/2))
    return wf_ini + filter

def lens_aberration(x,x_abcen,alpha,focus,ap):
    x_cen = (x_abcen - 0.5)*ap
    return np.exp(1.j*alpha * 1e9 * pow((x-x_cen)/focus,3))

def lens_xre(x,x_a,focus_x = 1.2e-3, x_abcen = 0.5, defocus=400e-6,ap_x = 40e-6,wl = 7.29e-11,alpha_x = -0.05):
    x_cen = (x_abcen - 0.5)*ap_x
    k = 2*np.pi / wl
    M_x = (defocus + focus_x) / defocus / focus_x
    z = defocus + focus_x
    f =lens_aberration(x=x_a,x_abcen=x_abcen,focus=focus_x,alpha=alpha_x,ap=ap_x) * np.exp(1.j* k/2/z*x_a**2) * np.exp(-1.j*k/M_x * x *x_a)

    return np.real(f)

def lens_xim(x,x_a,focus_x = 1.2e-3, x_abcen = 0.5, defocus=400e-6,ap_x = 40e-6,wl = 7.29e-11,alpha_x = -0.05):
    x_cen = (x_abcen - 0.5)*ap_x
    k = 2*np.pi / wl
    M_x = (defocus + focus_x) / defocus / focus_x
    z = defocus + focus_x
    f =lens_aberration(x=x_a,x_abcen=x_abcen,focus=focus_x,alpha=alpha_x,ap=ap_x) * np.exp(1.j* k/2/z*x_a**2) * np.exp(-1.j*k/M_x * x *x_a)
    return np.imag(f)

def lens_xc(x,x_arr,ap_x = 40e-06,focus_x = 1.2e-3, defocus = 400e-6,x_abcen=0.5, alpha_x = -0.05, wl = 7.29e-11):
    k = 2*np.pi / wl
    M_x = (defocus + focus_x) / defocus / focus_x
    z = defocus + focus_x
    A = - np.exp(1.j*k*z) / wl /defocus/focus_x * np.exp(1.j * k/2 * x**2 / M_x)

    fn_x = np.int(ap_x**2 / wl / (focus_x + defocus))

    wf_lens_x = []
    for i in range(len(x_arr)):
        Re = lens_xre(x,x_arr[i])
        Im = lens_xim(x,x_arr[i])
        Re_inte = integrate.quad(Re,-ap_x / 2,ap_x / 2,epsabs = 10-7,epsrel = 10-9,limit = 1000 * fn_x)
        Im_inte = integrate.quad(Im,-ap_x / 2,ap_x / 2,epsabs = 10-7,epsrel = 10-9,limit = 1000 * fn_x)
        wf_lens_x.append((Re_inte+1.j*Im_inte)*A)
    return wf_lens_x

def lens_wp_y(y,y_arr,ap_y = 40e-06,focus_y = 1.2e-3, defocus = 400e-6,y_abcen = 0.5 , alpha_y = -0.05, wl = 7.29e-11):
    y_cen = (y_abcen - 0.5)*ap_y
    k = 2*np.pi / wl
    M_y = (defocus + focus_y) / defocus / focus_y
    z = defocus + focus_y


    def func_re(y,y_arr,focus_y = focus_y, y_abcen = y_abcen, defocus = defocus):
        y_cen = (y_abcen - 0.5)*ap_y
        z = defocus + focus_y
        M_y=(defocus + focus_y) / defocus / focus_y
        f = np.exp(1.j*alpha_y*1e9*pow((y_arr-y_cen)/focus_y,3)) * np.exp(1.j* k/2/z*y_arr**2) * np.exp(-1.j*k/M_y * y *y_arr)
        return np.real(f)

    def func_im(y,y_arr,focus_y = focus_y, y_cen = y_cen, M_y=(defoucs + focus_y) / defoucs / focus_y, z=defoucs + focus_y):
        f = np.exp(1.j*alpha_y*pow((y_arr-y_cen)/focus_y,3)) * np.exp(1.j* k/2/z*x_arr**2) * np.exp(-1.j*k/M_y * y *y_arr)
        return np.im(f)

    A = - np.exp(1.j*k*z) / wl /defocus/focus_y * np.exp(1.j * k/2 * y**2 / M_y)
    fn_y = (ap_y**2 / wl / (focus_y + defocus))

    Re_func = func_re(y,y_arr)
    Re_inte = integrate.quad(Re_func,-ap_y / 2,ap_y / 2,epsabs = 10-7,epsrel = 10-9,limit = 1000 * fn)
    Im_func = func_im(x,x_arr)
    Im_inte = integrate.quad(Im_func,-ap_y / 2,ap_y / 2,epsabs = 10-7,epsrel = 10-9,limit = 1000 * fn_y)

    return (Re_inte+1.j*Im_inte)*A

def wavefront_defocus(x_arr,y_arr,wf_obj):
    for i in range(len(y_arr)):
        for j in range(len(x_arr)):
            wf_obj[i][j] = lens_wp_x(x_arr[j],x_arr)*lens_wp_y(y_arr[i],y_arr)
    return wf_obj

def propagate_fresnel(x,y,x_arr,wf_obj,y_arr,det_dist = 14e-03,wl=7.29e-11,defocus = 400e-06):
    k = np.pi * 2 / wl
    A = np.exp(1.j*k*det_dist)/1.j/wl/det_dist*np.exp(1.j*k/2/det_dist)*(x**2+y**2)
    M = (det_dist+defocus)/det_dist/defocus


    f = np.exp(-1.j*k/2/M*(x_arr**2 + y_arr**2)) * np.exp(-1.j*k/det_dist*(x*x_arr+y*y_arr))


    Re_inte = integrate.dblquad(np.real(f),-scipy.integrate.Inf,scipy.integrate.Inf)

    Im_inte = integrate.dblquad(np.imag(f),-scipy.integrate.Inf,scipy.integrate.Inf)

    return A*(Re_inte+1.j*Im_inte)
