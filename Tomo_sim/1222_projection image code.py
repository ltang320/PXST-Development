import numpy as np
import scipy.integrate as integrate
from scipy.signal import fftconvolve
import os,sys
import yaml
import matplotlib.pyplot as plt
import math

with open('configure.yml','r') as conf_para:
    conf_para = yaml.load(conf_para,Loader=yaml.FullLoader)
#print(conf_para)

def sample_coordinate(pixelsize_x = 55e-06,pixelsize_y = 55e-06,fs_size = 2000,ss_size = 2000,focus_x = 1.2e-3,focus_y = 1.0e-3,defocus = 400e-6, det_dist = 14e-03, ap_x = 40e-06, ap_y= 40e-6,wl = 7.29e-11):

    xx_span = fs_size * pixelsize_x
    yy_span = ss_size * pixelsize_y

    x_span_obj =  2*ap_x / focus_x * defocus
    y_span_obj =  2*ap_y / focus_y * defocus

    n_x = int(ap_x / wl * 2* ap_x/2/np.sqrt((ap_x/2)**2 + focus_x**2))*4
    n_y = int(ap_y / wl * 2* ap_y/2/np.sqrt((ap_y/2)**2 + focus_y**2))*4

    n_xx = int(x_span_obj/ wl * 2* ap_x/2/np.sqrt((ap_x/2)**2 + focus_x**2))*2
    n_yy = int(y_span_obj/ wl * 2* ap_y/2/np.sqrt((ap_y/2)**2 + focus_y**2))*2
    print(n_xx)
    nx_arr = np.arange(-n_x//2,n_x//2)
    ny_arr = np.arange(-n_y//2,n_y//2)
    ndecx_arr = np.arange(-n_xx//2, n_xx//2)
    ndecy_arr = np.arange(-n_yy//2, n_yy//2)

    delta_xlen = 2*ap_x / n_x
    delta_ylen = 2*ap_y / n_y
    #delta_xobj = ap_x/focus_x*defocus/n_x
    #delta_yobj = ap_y/focus_y*defocus/n_y
    delta_xobj = x_span_obj/n_x
    delta_yobj = y_span_obj/n_y
    delta_xdet = xx_span/n_xx
    delta_ydet = yy_span/n_yy

    return delta_xlen,delta_ylen,delta_xobj,delta_yobj,delta_xdet,delta_ydet,n_x,n_y,n_xx,n_yy,nx_arr,ny_arr,ndecx_arr,ndecy_arr

def wf_ini1d(corordinate,wavenumber,alpha,focus):
    u0 = np.exp(1.j*wavenumber/2/focus*corordinate**2 + 1e9j*alpha*(corordinate/focus)**3)
    return u0

def cov_kernel1d(corordinate,wavenumber,distance,af = 1,wl=7.29e-11):
    hn = distance/1.j/np.sqrt(wl)*np.exp(-1.j*wavenumber*np.sqrt(corordinate**2 + distance**2))/(corordinate**2 + distance**2)**0.75
    return hn

def wm(corordinate,n,af = 1):
    wm = np.exp(1.j*np.pi/n*corordinate**2*af)
    return wm

def hm(corordinate,n,af = 1):
    hm = np.exp(-1.j*np.pi/n*corordinate**2*af)
    return hm

alpha_x = conf_para['Lens']['alpha_x']
alpha_y = conf_para['Lens']['alpha_y']
focus_x = conf_para['Lens']['focus_x']
focus_y = conf_para['Lens']['focus_y']

det_dist = conf_para['Exp_geom']['det_dist']
defocus = conf_para['Exp_geom']['defocus']
ap_x = conf_para['Lens']['ap_x']
ap_y = conf_para['Lens']['ap_y']
wl = conf_para['Source']['wl']
k = 2*np.pi/wl

delta_xlen,delta_ylen,delta_xobj,delta_yobj,delta_xdet,delta_ydet,n_x,n_y,n_xx,n_yy,nx_arr,ny_arr,ndecx_arr,ndecy_arr = sample_coordinate(pixelsize_x = 55e-06,pixelsize_y = 55e-06,fs_size = 2000,ss_size = 2000,focus_x = 1.2e-3,focus_y = 1.0e-3,defocus = 400e-6, det_dist = 14e-03, ap_x = 40e-06, ap_y= 40e-6,wl = 7.29e-11)

#print(n_x)

x_len_arr = nx_arr * delta_xlen
u0_x = wf_ini1d(x_len_arr,wavenumber = k,alpha = alpha_x,focus=focus_x)
u0_x[np.abs(x_len_arr)>ap_x/2]=0
#%matplotlib widget
#plt.plot(x_len_arr,np.abs(u0_x))
#plt.plot(x_len_arr,np.unwrap(np.angle(u0_x)))
#print(n_x,u0_x)
af_obj_x = 1
af_len_x = defocus/focus_x
print(af_len_x)
x_obj_arr = nx_arr*delta_xobj

hn_x = cov_kernel1d(x_len_arr,wavenumber=k,wl=wl,distance = (focus_x+defocus),af = 1)
wm_x_len = wm(nx_arr,n=n_x,af=1)
hm_x_len = hm(nx_arr,n=n_x,af=1)
#%matplotlib widget
#plt.plot(wm_x_len)
wm_x_obj = wm(nx_arr,n=n_x,af=x_obj_arr)

u0_xf = np.fft.fftshift(np.fft.fft(u0_x))
hn_xf = np.fft.fftshift(np.fft.fft(hn_x))
#print(u0_xf.shape,hn_xf.shape,wm_x_len.shape)
gm_xf = u0_xf*hn_xf*wm_x_len
#print(hm_x_len.shape,gm_xf.shape)
u1_x = wm_x_obj * fftconvolve(gm_xf,hm_x_len,mode='same')
#u1_try = fftconvolve(u0_xf*hn_xf,hm_x_len,mode='same')
#%matplotlib widget
#plt.plot(np.abs(gm_xf))

y_len_arr = ny_arr * delta_ylen
u0_y = wf_ini1d(y_len_arr,wavenumber = k,alpha = alpha_y,focus=focus_y)
u0_y[np.abs(y_len_arr)>ap_y/2]=0
#af_obj_y = delta_ylen*delta_yobj*n_y/wl/(focus_y+defocus)
af_obj_y =1
af_len_y = defocus/focus_y
print(af_obj_y)
hn_y = cov_kernel1d(y_len_arr,wavenumber=k,wl=wl,distance = (focus_y+defocus),af = af_obj_y)

wm_y_len = wm(ny_arr,n=n_y,af=af_obj_y)
hm_y_len = hm(ny_arr,n=n_y,af=af_obj_y)

y_obj_arr = ny_arr*delta_yobj
wm_y_obj = wm(ny_arr,n=n_y,af=af_len_y)
u0_yf = np.fft.fftshift(np.fft.fft(u0_y))
hn_yf = np.fft.fftshift(np.fft.fft(hn_y))
gm_yf = u0_yf*hn_yf*wm_y_len

u1_y = wm_y_obj * fftconvolve(gm_yf,hm_y_len,mode='same')
#%matplotlib widget
#plt.plot(np.abs(gm_yf))
u1_x = np.reshape(u1_x,(len(u1_x),1))
u1_y = np.reshape(u1_y,(1,len(u1_y)))
u1_obj = u1_x * u1_y
