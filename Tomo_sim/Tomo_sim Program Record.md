footer: Program Record
slidenumbers: true
autoscale: true

# Tomography Simulation Project

- Build Class tomoSim

	- def __init__(self, tomo_sim_params, bsteps=none)

	- Here the bsteps is related to barcode. I should modify it.

	- def __getattr__(self,attr)

	- def _init_logging(self)

	- def _init)coord(self)

---

- Generate the initial coordinate arrays

	- ap_x and ap_y belongs to the lens parameters

	- xx_span = fs_size * pix_size (eg. fs_size = 2000)

	- yy_span = ss_size * pix_size (eg. xx_size = 1000)

	- x_span = 1.6*ap_x/focus/*defocus (According to similar triangle)

	- y_span = 1.2*ap_y

	- x_arr, y_arr, xx_arr, yy_arr: x axis and y axis value

- def_init_lens(self):

	- initialize wave fields at the sample's plane

	- self.logger.info generating wave fields at the sample's plane

---

- Term specify:

	- wl: wavelength

	- energy = e_to_wl = 1.2398419843320026e-06

	- det_dist sample to detector distance

---

# List of experimental parameters

## Experimental geometry parameters

- defocus  Lens' defocus distance [um]

- det_dist : Distance between the barcode and the detector [um]

- step_size : Scan step size [um]

- n_frames : Number of frames

---

## Detector parameters:

- fs_size : Detector's size along the fast axis in pixels

- ss_size : Detector's size along the slow axis in pixels

- pix_size : Detector's pixel size [um]

--- 

## Source parameters:

- p0 : Source beam flux [cnt/s]

- wl : Incoming beam's wavelength [um]

- th_s : Source rocking curve width [rad]

---

## Lens parameters:

- ap_x : Lens' aperture size along the x axis [um]

- ap_y : Lens' aperture size along the y axis [um]

- focus : Focal distance [um]

- alpha : Third order abberations ceofficient [rad/mrad^3]

- x0 : Lens' abberations center point [0.0 - 1.0]

--- 

## Barcode sample parameters:

- bar_size : Average bar's size [um]

- bar_sigma : Bar bluriness width [um]

- bar_atn : Bar's attenuation coefficient [0.0 - 1.0]

- bulk_atn : Barcode's bulk attenuation coefficient [0.0 - 1.0]

- rnd_dev : Bar's coordinates random deviation [0.0 - 1.0]
	
- offset : Barcode's offset at the beginning and at the end of the scan from the detector's bounds [um]



---

# Lens part 

- lens_re and lens_im is the same
- In function lens_wp_c, re and im part integral and then multiply the cos(ph)+1j*sin(ph)

```python
cdef double lens_re(double xx, void* params) nogil:
    cdef:
        double x = (<double*> params)[0], wl = (<double*> params)[1]
        double f = (<double*> params)[2], df = (<double*> params)[3]
        double a = (<double*> params)[4], xc = (<double*> params)[5]
        double ph, ph_ab
    ph = -pi * xx**2 / wl * df / f / (f + df) - 2 * pi / wl / (f + df) * x * xx
    ph_ab = -a * 1e9 * ((xx - xc) / f)**3
    return cos(ph + ph_ab)

cdef double lens_im(double xx, void* params) nogil:
    cdef:
        double x = (<double*> params)[0], wl = (<double*> params)[1]
        double f = (<double*> params)[2], df = (<double*> params)[3]
        double a = (<double*> params)[4], xc = (<double*> params)[5]
        double ph, ph_ab
    ph = -pi * xx**2 / wl * df / f / (f + df) - 2 * pi / wl / (f + df) * x * xx
    ph_ab = -a * 1e9 * ((xx - xc) / f)**3
    return sin(ph + ph_ab)

cdef complex_t lens_wp_c(double x, double wl, double ap, double f,
                       double df, double a, double xc) nogil:
    cdef:
        double re, im, ph = pi / wl / (f + df) * x**2
        double params[6]
        int fn = <int> (ap**2 / wl / (f + df))
        gsl_function func
    params[0] = x; params[1] = wl; params[2] = f
    params[3] = df; params[4] = a; params[5] = xc
    func.function = &lens_re; func.params = params
    re = gsl_quad(func, -ap / 2, ap / 2, 1e-9, 1e-7, 1000 * fn)
    func.function = &lens_im
    im = gsl_quad(func, -ap / 2, ap / 2, 1e-9, 1e-7, 1000 * fn)
    return (re + 1j * im) * (cos(ph) + 1j * sin(ph))
```

---

# Two methods of propagation in the objects

- The exit-surface at the lens' plane:

$$
U_0(x_0) = \Pi (a_x x_0) \exp
        \left[ -\frac{j \pi x_0^2}{\lambda f} + j \alpha
        \left( \frac{x_0 - x_c}{f} \right)^3 \right]
$$

---

- wavefront : $$U_0$$ propagates to the sample's plane which
    is $$f + z_1$$ downstream from the lens. According to
    the Fresnel diffraction theory (without the normalizing
    coefficient before the integral):

$$
U(x) = \int_{-a_x / 2}^{a_x / 2}
        e^{-\frac{j k z_1 x_0^2 }{2f(z_1 + f)}}
        e^{j\alpha\left(\frac{x_0 - x_c}{f}\right)^3} 
        e^{j\frac{2 \pi}{\lambda z} x x_0} dx_0
$$ 

---

# The Fraunhofer integral transform is defined as (without the normalizing coefficient before the integral)

$$
U(x) = e^{-\frac{j k x^2}{2 z}} \int_{-\infty}^{+\infty}
        U_0(x_0) e^{j\frac{2 \pi}{\lambda z} x x_0} dx_0
$$

---

# The exit-surface at the lens' plane:

$$
U_0(x_0,y_0) = \Pi (a_x x_0,a_y y_0)\exp[-\frac{j \pi x_0^2}{\lamda f_1} + j \alpha \frac{x_0 - x_c}{f}^3]\e1xp[-\frac{j \pi y_0^2}{\lamda f_2} + j \alpha \frac{y_0 - y_c}{f}^3]

$$

---

# You can also turn the filter off.

![left filtered 250%](http://deckset-assets.s3-website-us-east-1.amazonaws.com/colnago2.jpg)


