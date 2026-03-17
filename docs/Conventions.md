# Conventions
We define the mass ratio as $q \equiv m_1/m_2 \ge 1$, where $m_{1,2}$ are the masses of the two black holes, the reduced mass of the system as $\mu = m_1 m_2 /M$, with $M = m_1 + m_2$ being the total mass of the system, and the symmetric mass ratio as $\nu \equiv m_1 m_2/M^2 = \mu/M$. 

We evolve dimensionless Cartesian coordinates $\vec{r} \equiv \vec{R}/M = \\{ x, y, z \\}$ and momenta $\vec{p} \equiv \vec{P}/\mu = \\{ p_x, p_y, p_z \\}$ with dimensionless time $t \equiv T/M$, in geometric units $G = c = 1$. We further set $M = 1$.

The individual dimensionless spins of the black holes are defined as $\vec{\chi}\_i \equiv \vec{S}\_i / m_i^2$ for $i = 1, 2$, and $\vec{\chi}\_i = \\{ \chi_{i,x}, \chi_{i,y}, \chi_{i,z} \\}$.  

## Input parfile

**Binary parameters**

| Parameter | Description |
| --------- | ----------- |
| `q`       | Mass ratio |
| `chi1x0`  | Initial $x$ component of the (dimensionless) spin of the first black hole |
| `chi1y0`  | Initial $y$ component of the (dimensionless) spin of the first black hole |
| `chi1z0`  | Initial $z$ component of the (dimensionless) spin of the first black hole |
| `chi2x0`  | Initial $x$ component of the (dimensionless) spin of the second black hole |
| `chi2y0`  | Initial $y$ component of the (dimensionless) spin of the second black hole |
| `chi2z0`  | Initial $z$ component of the (dimensionless) spin of the second black hole |

**Orbital parameters**

| Parameter | Description |
| --------- | ----------- |
| `x0`      | $x$ component of the initial separation of the binary |
| `y0`      | $y$ component of the initial separation of the binary |
| `z0`      | $z$ component of the initial separation of the binary |

**Constants of motion**

| Parameter | Description |
| --------- | ----------- |
| `E`       | Energy |
| `lz`      | $z$ component of the (dimensionless) orbital angular momentum, $\vec{l} = \vec{r} \times \vec{p}$ |

**Options for different functions**

| Parameter | Description | Options |
| --------- | ----------- | ------- |
| `pots`    | Choice for the orbital potentials | `"resummed"`, `"non-resummed"` |

**ODE solver settings**

| Parameter | Description | Options |
| --------- | ----------- | ------- |
| `dt`      | Time step (when there is no time transformation - see Model.md) | - |
| `tmax`    | Maximum integration time | - |
| `solver`  | Type of ODE solver | `"RK4"`, `"RK-GL6"`|
| `step`    | Type of time step (see Model.md) | `"standard"`, `"transformed"`|

## Output

The code outputs two files, `dyn.txt` and `metadata.txt`. The latter specifies the given initial conditions, together with other choices for the dynamics and for the integration scheme. The first file has the following columns:

| # | Column    | Description | 
| - | --------- | ----------- |
| 1 | `t`       | time |
| 2 | `x`       | $x$ component of the separation of the binary |
| 3 | `y`       | $y$ component of the separation of the binary |
| 4 | `z`       | $z$ component of the separation of the binary |
| 5 | `px`      | $x$ component of the conjugate momentum $\vec{p}$ |
| 6 | `py`      | $y$ component of the conjugate momentum $\vec{p}$ |
| 7 | `pz`      | $z$ component of the conjugate momentum $\vec{p}$ |
| 8 | `chi1x`   | $x$ component of the (dimensionless) spin of the first black hole |
| 9 | `chi1y`   | $y$ component of the (dimensionless) spin of the first black hole |
| 10 | `chi1z`  | $z$ component of the (dimensionless) spin of the first black hole |
| 11 | `chi2x`  | $x$ component of the (dimensionless) spin of the second black hole |
| 12 | `chi2y`  | $y$ component of the (dimensionless) spin of the second black hole |
| 13 | `chi2z`  | $z$ component of the (dimensionless) spin of the second black hole |
| 14 | `E`      | Real energy of the binary |
| 15 | `lx`     | $x$ component of the (dimensionless) orbital angular momentum ($\vec{l} = \vec{r} \times \vec{p}$) |
| 16 | `ly`     | $y$ component of the (dimensionless) orbital angular momentum ($\vec{l} = \vec{r} \times \vec{p}$) |
| 17 | `lz`     | $z$ component of the (dimensionless) orbital angular momentum ($\vec{l} = \vec{r} \times \vec{p}$) |

