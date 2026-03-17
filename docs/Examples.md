## Examples

Below we showcase some results that can be obtained with the code. In particular, we make use of the parameter files that can be found in `C/parfiles`. For each configuration, we specify the starting values of the binary parameters (mass ratio and spins), the desired type of initial conditions, and the orbital parameters. There are three types of initial conditions:
- Circular initial conditions: with given initial radius along the $x$ axis ($x_0$), the code finds the initial angular momentum along the $z$ direction, $p_\varphi = l_z$, as the root of $\partial H_{\rm eff} / \partial x$, and sets $p_{y,0} = p_{\varphi,0}/x_0$. 
- Equatorial eccentric initial conditions: for given inital energy and angular momentum along the $z$ direction ($E_0, p_{\varphi,0}$), the code finds the equatorial turning points and starts the motion at the apoapsis.
- Generic initial conditions: for given initial orbital parameters and inital energy and angular momentum along the $z$ direction, the code sets $p_{y,0} = p_{\varphi,0}/x_0$ and finds $p_{z,0}$ as the value that yields $H_{\rm eff} = E_0$.

All the above initial conditions exploit the effective Hamiltonian $H_{\rm eff}$ as an effective potential. In all the three cases, one can set the spins either only along the $z$ axis, or with in-plane components as well. Therefore, we can have six types of orbits, displayed below. For the orbits for which it is relevant, we also plot the evolution of the (dimensionless) angular momentum $\vec{l} = \vec{r} \times \vec{p}$ and of the (dimensionless) spins, $\vec{\chi}_1, \vec{\chi}_2$, multiplied by 3 for visual purposes.

### Circular initial conditions, aligned spins

```
# binary parameters
q      = 10.  
chi1x0 = 0.
chi1y0 = 0. 
chi1z0 = 0.65
chi2x0 = 0.
chi2y0 = 0. 
chi2z0 = 0.65

# type of initial conditions
motion = circular

# orbital parameters
x0 = 15.
y0 = 0.
z0 = 0.
```

| Trajectory | Energy conservation |
| :---: | :---: |
| <img src="figs/circ_q_10.00_chi1z0_0.65_chi2z0_0.65_traj.png" height="400"> | <img src="figs/circ_q_10.00_chi1z0_0.65_chi2z0_0.65_energycons.png" height="400"> |

### Circular initial conditions, misaligned spins

```
# binary parameters
q      = 10.
chi1x0 = 0.16
chi1y0 = 0.29
chi1z0 = 0.56
chi2x0 = 0.47
chi2y0 = 0.51
chi2z0 = 0.65

# type of initial conditions
motion = circular

# orbital parameters
x0 = 15.
y0 = 0.
z0 = 0.
```

| Trajectory | Top view |
| :---: | :---: |
| <img src="figs/circ_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_traj.png" height="400"> | <img src="figs/circ_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_traj_topview.png" height="400"> |

| Energy conservation | Angular momentum and spins evolution |
| :---: | :---: |
| <img src="figs/circ_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_energycons.png" height="400"> | <img src="figs/circ_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_angmomandspins.png" height="400"> |

### Equatorial initial conditions, aligned spins

```
# binary parameters
q      = 10.  
chi1x0 = 0.
chi1y0 = 0.
chi1z0 = 0.65
chi2x0 = 0.
chi2y0 = 0.
chi2z0 = 0.65

# type of initial conditions
motion = equatorial-eccentric

# orbital parameters
E0  = 0.965
pphi0 = 3.8
```

| Trajectory | Energy conservation |
| :---: | :---: |
| <img src="figs/eq_q_10.00_chi1z0_0.65_chi2z0_0.65_traj.png" height="400"> | <img src="figs/eq_q_10.00_chi1z0_0.65_chi2z0_0.65_energycons.png" height="400"> |

### Equatorial initial conditions, misaligned spins

```
# binary parameters
q      = 10.  
chi1x0 = 0.16
chi1y0 = 0.29
chi1z0 = 0.56
chi2x0 = 0.47
chi2y0 = 0.51
chi2z0 = 0.65

# type of initial conditions
motion = equatorial-eccentric

# orbital parameters
E0  = 0.97
pphi0 = 3.8
```

| Trajectory | Top view |
| :---: | :---: |
| <img src="figs/eq_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_traj.png" height="400"> | <img src="figs/eq_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_traj_topview.png" height="400"> |

| Energy conservation | Angular momentum and spins evolution |
| :---: | :---: |
| <img src="figs/eq_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_energycons.png" height="400"> | <img src="figs/eq_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_angmomandspins.png" height="400"> |

### Generic initial conditions, aligned spins

```
# binary parameters
q      = 10.  
chi1x0 = 0.
chi1y0 = 0.
chi1z0 = 0.65
chi2x0 = 0.
chi2y0 = 0.
chi2z0 = 0.65

# type of initial conditions
motion = generic

# orbital parameters
x0 = 10. 
y0 = 0.
z0 = 0.
E0    = 0.97
pphi0 = 3.
```

| Trajectory | Energy conservation |
| :---: | :---: |
| <img src="figs/gen_q_10.00_chi1z0_0.65_chi2z0_0.65_traj.png" height="400"> | <img src="figs/gen_q_10.00_chi1z0_0.65_chi2z0_0.65_energycons.png" height="400"> |

| Angular momentum and spins evolution |
| :---: |
| <img src="figs/gen_q_10.00_chi1z0_0.65_chi2z0_0.65_angmomandspins.png" height="400"> |

### Generic initial conditions, misaligned spins

```
# binary parameters
q      = 10.  
chi1x0 = 0.16
chi1y0 = 0.29
chi1z0 = 0.56
chi2x0 = 0.47
chi2y0 = 0.51
chi2z0 = 0.65

# type of initial conditions
motion = generic

# orbital parameters
x0 = 10.
y0 = 0.
z0 = 0.
E0    = 0.97
pphi0 = 3.5
```
| Trajectory | Energy conservation |
| :---: | :---: |
| <img src="figs/gen_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_traj.png" height="400"> | <img src="figs/gen_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_energycons.png" height="400"> |

| Angular momentum and spins evolution |
| :---: |
| <img src="figs/gen_q_10.00_chi1-0_0.16_0.29_0.56_chi2-0_0.47_0.51_0.65_angmomandspins.png" height="400"> |




