# cubrium
A **cub**e in equilib**rium**

[![PyPI version shields.io](https://img.shields.io/pypi/v/cubrium.svg)](https://pypi.python.org/pypi/cubrium/)
![Code coverage](coverage.svg)
![Made with love in Graz](https://madewithlove.now.sh/at?heart=true&colorA=%233b3b3b&colorB=%231f744f&text=Graz)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

<img src="https://raw.githubusercontent.com/adtzlr/cubrium/main/cube.png" width="75%">

Cubrium is a toolbox for the definition, solution and post-processing of homogenous loadcases in continuum mechanics of solids (statics). It uses [contique](https://github.com/adtzlr/contique/) for the numeric continuation of the nonlinear equilibrium equations. If you use `cubrium` or `contique` in your scientific publications, please cite my work. I'll provide a citable template with a DOI in the future.

## Example 101 a.k.a `hello cubrium` 😎
This is an example which solves a cube with a Saint Venant-Kirchhoff (SVK) material for the case of uniaxial loading. In the first step, we `init` a model.

```python
import cubrium

MDL = cubrium.init()
```

### Constitution
In the second part, we have to define the constitutive law. We can either use one of `cubrium`'s models or use our own `umat` (user material). This time we're using our own `umat` function for a simple SVK material.

<img src="https://latex.codecogs.com/gif.latex?\boldsymbol{S}&space;=&space;2&space;\mu&space;\&space;\boldsymbol{E}&space;&plus;&space;\gamma&space;\&space;\mathrm{tr}(\boldsymbol{E})&space;\&space;\boldsymbol{1}" title="\boldsymbol{S} = 2 \mu \ \boldsymbol{E} + \gamma \ \mathrm{tr}(\boldsymbol{E}) \ \boldsymbol{1}" />

```python
import numpy as np

def umat_svk(F, parameters):
    """(U)ser (MAT)erial Function.
    Returns First Piola-Kirchhoff stress tensor for a given
    deformation gradient tensor with a list of material parameters."""

    # expand list of material parameters
    mu, K = parameters[:2]
    gamma = K - 2 / 3 * mu

    C = F.T @ F
    E = 1 / 2 * (C - np.eye(3))
    
    S = 2 * mu * E + gamma * np.trace(E) * np.eye(3)

    return F @ S
```

Now we have to link our `umat` to the `cubrium` model definition and specify material parameters.

```python
MDL.GLO.constitution.umat = umat_svk
MDL.GLO.constitution.parameters = [1.0, 5000.0]
```

### Loadcase (Kinematics and Kinetics)
A loadcase is defined with exactly **9** equations for the unsymmetric or **6** equations for the full-symmetric case. This contains either kinematic or kinetic types of equations. For the case of uniaxial loading we are building this loadcase for ourselfes. We apply an external normal force 1 and set all external shear forces and normal forces 2 and 3 to zero. A symmetric solution is enforced (no rigid body rotation is allowed). The load-proportionaly-factor is applied to the normal forces (`lpftype=0`). Finally we specify a `title` for the loadcase. This will later effect the output filenames.

```python
def uniaxial(MDL):
    MDL.EXT.force.normal[0] = 1
    MDL.EXT.force.normal[1] = 0
    MDL.EXT.force.normal[2] = 0

    MDL.EXT.force.shear[0, 1] = 0
    MDL.EXT.force.shear[1, 2] = 0
    MDL.EXT.force.shear[0, 2] = 0

    MDL.EXT.gridvec.symmetry = [1, 1, 1]

    MDL.GLO.lpftype = 0
    MDL.GLO.title = "Uniaxial"
    return MDL
```

Again, we have to link our loadcase to the `cubrium` model and update the model with the new loadcase settings.

```python
MDL = uniaxial(MDL)
MDL = cubrium.update(MDL)
```

### Solver
Starting from a valid initial solution everything is ready to solve the model. **Hint**: `x0` are the flattened components of the displacement gradient w.r.t. the undeformed coordinates (=primary unknows of the problem).

```python
Res = cubrium.solve(MDL)(
    x0   = np.zeros(9),
    lpf0 = 0.0,
)
```

The results contain the extended unknowns `y = (x, lpf)` but no information about the internal quantities of the model. Therefore we extract the extended unknows from the Result object (`Res`) and recover these internal quantities (e.g. reaction forces) for all steps.

```python
Y = np.array([res.x for res in Res])
history = cubrium.recover(Y, MDL)
```

### Plots and Post-processing
We plot the axial stretch vs. load-proportionality-factor in direction 1.

```python
import matplotlib.pyplot as plt

plt.plot(1+Y[:, 0], Y[:, -1], "-")
plt.xlabel("stretch $\lambda_1$")
plt.ylabel("load-proportionality-factor LPF")
```

<img src="https://raw.githubusercontent.com/adtzlr/cubrium/main/scripts/Uniaxial_stretch-lpf.svg" width="75%">

Using `meshio` we are able to export our solution in the `xdmf` file format which may be further post-processed by ParaView.

```python
cubrium.writer.xdmf(
        history,
        filename = MDL.GLO.title,
    )
```

An exemplary scene for ParaView 5.9.0 is available to [download](https://raw.githubusercontent.com/adtzlr/cubrium/main/scripts/paraviewstatecube.pvsm). Import it in ParaView (File - Load state) and choose "Choose File Names" as shown below. Voilà, a nice cube scene in 3D with a cube colored in "Cauchy Stress XX" and reaction forces scaled and colored in "Reaction Force Magnitude" is ready to animate. The whole script of this example may be [downloaded here](https://raw.githubusercontent.com/adtzlr/cubrium/main/scripts/script101_hellocubrium.py).

<img src="https://raw.githubusercontent.com/adtzlr/cubrium/main/scripts/paraviewimportscene.png" width="75%">

<img src="https://raw.githubusercontent.com/adtzlr/cubrium/main/scripts/paraviewcube.png" width="75%">

Finally you can watch the animated cube during the deformation.

<a href="ttps://raw.githubusercontent.com/adtzlr/cubrium/main/scripts/script101_hellocubrium_video.ogv"><img src="https://raw.githubusercontent.com/adtzlr/cubrium/main/scripts/script101_hellocubrium_video.gif" href="" width="75%"></a>

Have fun using `cubrium`! If you find any bugs please submit an issue.