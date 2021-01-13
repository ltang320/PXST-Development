import os,sys
import numpy as np
import yaml
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math

"""
------------
Parameters

"""
with open('configure.yml','r') as conf_para:
    conf_para = yaml.load(conf_para,Loader=yaml.FullLoader)

"""
------------
wavefront_initialize

>> input
pixelsize_x = 55e-06,pixelsize_y=55e-06,
fs_size = 2000,ss_size = 2000,
focus_x = 1.2e-3,focus_y = 1.0e-3,
defocus = 400e-6,
det_dist = 14e-03,
ap_x = 40e-06, ap_y= 40e-6,
wl = 7.29e-11,
amplitude_value=0.0

>> output
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

    wf_obj = np.zeros((n_x,n_y),dtype='complex')

    # coordinate in detector plan
    xx_arr = np.linspace(-xx_span / 2, xx_span / 2, fs_size, endpoint=False)
    yy_arr = np.linspace(-yy_span / 2, yy_span / 2, ss_size, endpoint=False)

    return x_arr,y_arr,xx_arr,yy_arr,wf_obj

"""
lens wavefront

Parameters:
------------
r : coordinates
f : focus of lens
df: defocus of the object
a : alpha, Third order abberations coefficient [rad/mrad^3]
cen_ab : center point of the lens' abberations

output
------
wavefront_lens,err

"""
def lens_wf(x_arr, y_arr, wf_obj, ap_x = 40e-06,ap_y = 40e-06, focus_x = 1.2e-3, focus_y=1.0e-3, x_abcen = 0.5, y_abcen = 0.5, alpha_x = -0.05, alpha_y = -0.05, wl = 7.29e-11,defocus =400e-06):
    xx = x_arr.copy()
    yy = y_arr.copy()
    
    wf_lens = np.array(np.meshgrid(y_arr,x_arr))
    wf_obj_cor = np.array(np.meshgrid(yy,xx))
    wavefront_lens = np.zeros_like(wf_obj,dtype='complex')

    wavenumber = 2*np.pi / wl
    
    z_dis = focus_y + defocus

    M_x = (focus_x+defocus)/focus_x
    M_y = (focus_y+defocus)/focus_y

    A = wavenumber/1.j/2/np.pi/z_dis

    ph_0 = wavenumber* 1.j / 2 / z_dis * (wf_obj_cor[0,:,:]**2 + wf_obj_cor[1,:,:]**2) + 1.j*wavenumber*z_dis

    x_cen = (x_abcen - 0.5)*ap_x
    y_cen = (y_abcen - 0.5)*ap_y

    ph_x = -wavenumber / 2 / M_x / focus_x * wf_lens[0,:,:]**2
    ph_ab_x = alpha_x * 1e9 * ((wf_lens[0,:,:] - x_cen) / focus_x) **3

    ph_y = -wavenumber / 2 / M_y / wf_lens[1,:,:] * y_arr**2
    ph_ab_y= alpha_y * 1e9 * ((wf_lens[1,:,:] - y_cen) / focus_y) **3

    ph_mix = wavenumber / defocus * (wf_obj_cor[0,:,:]*wf_lens[0,:,:] + wf_obj_cor[1,:,:]*wf_lens[1,:,:])

    func = np.exp(1.j * (ph_x + ph_ab_x + ph_y + ph_ab_y + ph_mix))

    for i in range(y_arr.size):
        for j in range(x_arr.size):
            wavefront_lens[i][j], err = integrate.dblquad(func,
                                            -ap_x / 2,
                                            ap_x / 2,
                                            -ap_y / 2,
                                            ap_y / 2,
                                            args=(),
                                            epsabs=1e-07,
                                            epsrel=1e-09)

    wavefront_lens *= A*np.exp(ph_0)
    
    return wavefront_lens,err

def propagator2d_integrate(x_arr,y_arr,xx_arr,yy_arr,wf_dec,wavefront_lens, det_dist = 14e-03, wl = 7.29e-11 ):

    # convolving with the Fresnel kernel via FFT multiplication
    p_xy = np.array(np.meshgrid(y_arr,x_arr))
    det_xy = np.array(np.meshgrid(yy_arr,xx_arr))
    #wf_propagated = wavefront_lens
    wf_progagated = np.zeros_like(wf_dec,dtype='complex')

    wavenumber = 2 * np.pi / wl

    ph = wavenumber / 2 / det_dist

    for i in range(yy_arr.size):
        for j in range(xx_arr.size):
            ph_x = wavenumber/ det_dist * p_xy[0,:,:] * det_xy[0,j,i]
            ph_y = wavenumber/ det_dist * p_xy[1,:,:] * det_xy[0,j,i]
            value = wavefront_lens * np.exp(-ph_x-ph_y)
            wf_propagated[i][j] *= np.exp(1.j*ph) * integrate.simps(integrate.simps(value,ph_y),ph_x)

    return wf_propagated

x_arr,y_arr,xx_arr,yy_arr,wf_dec = wavefront_initialize(pixelsize_x = 55e-06,pixelsize_y = 55e-06,fs_size = 20,ss_size = 20,focus_x = 1.2e-3,focus_y = 1.0e-3,defocus = 400e-6, det_dist = 14e-03, ap_x = 40e-06, ap_y= 40e-6,wl = 7.29e-4,amplitude_value=0.0)

wavefront_lens,err = lens_wf(x_arr, y_arr, xx_arr, yy_arr, wf_dec)

wf_progagated = propagator2d_integrate(x_arr,y_arr,xx_arr,yy_arr,wf_dec,wavefront_lens, det_dist = 14e-03, wl = 7.29e-11 )

fig,(ax1, ax2) = plt.subplots(1,2)
ax1.set_title('amplitude')
im1 = ax1.imshow(np.real(wf_progagated))

ax2.set_title('phase')
im2 = ax2.imshow(np.imag(np.unwrap(wf_progagated)))

plt.tight_layout()
# Make space for title
plt.subplots_adjust(top=0.85)
plt.show()
