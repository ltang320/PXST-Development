footer: Weekly Report
slidenumbers: true
autoscale: true

# Week 50 1207 - 1213 Review

- Write and adjust the program of cone-beam 2D propagation code
- Consider the procedure to generate the 3D projection image


---

# Core procedure of propagation -1 

- Generate the initial coordinate and wavefront

- Generate lens_wavefront
	- Considered that the variables in x-axis and y-axis could seperate

	$$
	U_0( \zeta , \eta) = \exp[-\frac{ik}{2f_x}x_0^2 + i \alpha_x (\frac{x_0 - x_c}{f_x})^3]  \exp [-\frac{ik}{2f_y}y_0^2 + i \alpha_y (\frac{y_0 - y_c}{f_y})^3]
	$$

- Propagate the lens_wavefront to object surface

	- Fresnel diffraction to object surface 
	- Consider the aperture 

	$$
	U_{obj}(x,y) = \frac {\exp [jkz_{def}]}{i \lambda z_{def}} \exp [\frac{ik(x^2 + y^2)}{2(f_y+Z_{def})}y_0^2 \int_{\frac {-a_y}{2}}^{\frac {a_y}{2}} \int_{\frac {-a_x}{2}}^{\frac {a_x}{2}} \exp[-\frac{ik z_{def}}{2(f_x+Z_{def})}x_0^2 + i \alpha_x (\frac{x_0 - x_c}{f_x})^3]  \exp [-\frac{ik z_{def}}{2(f_y+Z_{def})}y_0^2 + i \alpha_y (\frac{y_0 - y_c}{f_y})^3] \exp [j \frac{k}{z_{def}}(x \zeta + y \eta)]d\zeta d\eta
	$$

---

# Core procedure of propagation -2

- Multiply the matrix with the object 

- Propagate to the detector plane

	$$
	U_det(x,y) = \frac {\exp [jkz]}{i \lambda z} \int \int T( \zeta , \eta) U_{obj}(\zeta , \eta)\exp[-\frac{ik 
[(x - \zeta)^2 + (y - \eta)^2]}{2z}d\zeta d\eta
	$$


- The relationship between the coordinate matrix and the detector matrix is a little confused which generates some problems in the code. 

---
 
# Next step

- Debug the 2D propagation code and generate the proper image
- Evaluate the calculation time of the code
- Consider the coordinate generation in 3D object
- Generate the projection image in one direction
- Rotate the object in different direction and generate the dataset
- Scanning different part of the object

