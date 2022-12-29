# Introduction to MaxwellWave

```@contents
Pages = [
    "index.md",
    "wave2D.md",
    "wave3D.md",
    "wave3DnonlinearChi3.md",
    "auxiliary.md",
]
```
# Maxwell Wave
[![Build Status](https://github.com/lheimabch/MaederHeimbach/actions/workflows/CI.yml/badge.svg)](https://github.com/lheimabch/MaederHeimbach/actions/workflows/CI.yml)

The goal of this project is to solve the electromagnetic wave equation in 3D.

The documentation is hosted under [MaxwellWave](https://lheimabch.github.io/MaederHeimbach/dev/)

## Authors
- Lothar Heimbach, Master's Student in Computational Science and Engineering at ETH Zurich
- Alexander Maeder, Master's Student in Electrical Engineering and Information Technology at ETH Zurich

# Derivation Partial Differentaial Equation
We start with the general Maxwell's equation:

$$\nabla \cdot \vec{D}\left(\vec{r},t\right) = p_0\left(\vec{r},t\right)$$

$$\nabla \times \vec{E}\left(\vec{r},t\right) = - \frac{\partial}{\partial t}\vec{B}\left(\vec{r},t\right)$$

$$\nabla \times \vec{H}\left(\vec{r},t\right) = \frac{\partial}{\partial t} \vec{D}\left(\vec{r},t\right) + \vec{j}_0\left(\vec{r},t\right)$$

$$\nabla \cdot \vec{B}\left(\vec{r},t\right) = 0$$

We have to further supplement this system of equation with constitutive relations.

## Linear Lossy Problem
We choose for the first problem setting a nonmagnetic dielectric with no free charges and no nonlinearities:<br />

$$p_0\left(\vec{r},t\right) = 0$$

$$\vec{j}_0\left(\vec{r},t\right) = \sigma\left(\vec{r},t\right) \vec{E}\left(\vec{r},t\right)$$

$$\vec{D}\left(\vec{r},t\right) = \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \vec{E}\left(\vec{r},t\right)$$

$$\vec{B}\left(\vec{r},t\right) = \mu_0 \vec{H}\left(\vec{r},t\right)$$

We then insert these constitutive relations in the maxwell equations and arrive at the system:

$$\nabla \cdot \left( \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \vec{E}\left(\vec{r},t\right) \right) = 0$$

$$\nabla \times \vec{E}\left(\vec{r},t\right) = - \mu_0 \frac{\partial}{\partial t}\vec{H}\left(\vec{r},t\right)$$

$$\nabla \times \vec{H}\left(\vec{r},t\right) = \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \frac{\partial}{\partial t} \vec{E}\left(\vec{r},t\right) + \sigma\left(\vec{r},t\right) \vec{E}\left(\vec{r},t\right)$$

$$\nabla \cdot \vec{H}\left(\vec{r},t\right) = 0$$

To further simplify the system, we apply multiple transformations and vector caculus identities:

$$\nabla \times \nabla \times \vec{E}\left(\vec{r},t\right) = - \mu_0 \frac{\partial}{\partial t} \nabla \times \vec{H}\left(\vec{r},t\right)$$

$$\nabla \times \nabla \times \vec{E}\left(\vec{r},t\right) = - \mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \frac{\partial^2}{\partial t^2}\vec{E}\left(\vec{r},t\right) -\mu_0 \left(\frac{\partial}{\partial t}\sigma\left(\vec{r},t\right)\right)\vec{E}\left(\vec{r},t\right) -\mu_0\sigma\left(\vec{r},t\right) \left(\frac{\partial}{\partial t}\vec{E}\left(\vec{r},t\right)\right)$$

$$\nabla \times \nabla \times \vec{E}\left(\vec{r},t\right) = \nabla \left(\nabla \cdot \vec{E}\left(\vec{r},t\right)\right) - \nabla^2 \vec{E}\left(\vec{r},t\right)$$

$$0 = \nabla \cdot \left(\varepsilon_r\left(\vec{r}\right) \vec{E}\left(\vec{r},t\right) \right) = \varepsilon_r\left(\vec{r}\right) \nabla \cdot \vec{E}\left(\vec{r},t\right) + \nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)$$

$$\nabla \cdot \vec{E}\left(\vec{r},t\right) = -\frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)}$$

$$\nabla \times \nabla \times \vec{E}\left(\vec{r},t\right) = -\nabla \left(\frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} \right) - \nabla^2 \vec{E}\left(\vec{r},t\right)$$

With all of these simplifications, we arrive at the final equation for the electric field, which we want to solve: 

$$\nabla^2 \vec{E}\left(\vec{r},t\right) + \nabla \left(\frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} \right) - \mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \frac{\partial^2}{\partial t^2}\vec{E}\left(\vec{r},t\right) -\mu_0 \left(\frac{\partial}{\partial t}\sigma\left(\vec{r},t\right)\right)\vec{E}\left(\vec{r},t\right) -\mu_0\sigma\left(\vec{r},t\right) \left(\frac{\partial}{\partial t}\vec{E}\left(\vec{r},t\right)\right) = 0$$

which can be equivalently formulated as: 

$$\frac{\partial}{\partial t}\begin{bmatrix}
        \vec{u}\left(\vec{r},t\right)\\
        \vec{v}\left(\vec{r},t\right)\\
    \end{bmatrix} 
    = \begin{bmatrix}
        0 & 1\\
        -\frac{\frac{\partial}{\partial t}\sigma\left(\vec{r},t\right)}{\varepsilon_0 \varepsilon_r\left(\vec{r}\right)} & -\frac{\sigma\left(\vec{r},t\right)}{\varepsilon_0 \varepsilon_r\left(\vec{r}\right)} \\
    \end{bmatrix}
    \begin{bmatrix}
        \vec{u}\left(\vec{r},t\right)\\
        \vec{v}\left(\vec{r},t\right)\\
    \end{bmatrix} + 
    \begin{bmatrix}
        0\\
        \frac{1}{\mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right)}\nabla^2 \vec{u}\left(\vec{r},t\right) + \frac{1}{\mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right)} \nabla \left(\frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{u}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} \right)\\
    \end{bmatrix}$$

$$\vec{u}\left(\vec{r},t\right) = \vec{E}\left(\vec{r},t\right)$$

$$\vec{v}\left(\vec{r},t\right) = \frac{\partial}{\partial t} \vec{E}\left(\vec{r},t\right)$$

## $\chi^3$ Nonlinear Lossy Problem
We choose for the second problem setting a more general setting such that third harmonic effects can be captured by the simulation, for this we assumse these constitutive relations:<br />

$$p_0\left(\vec{r},t\right) = 0$$

$$\vec{j}_0\left(\vec{r},t\right) = \sigma\left(\vec{r},t\right) \vec{E}\left(\vec{r},t\right)$$

$$\vec{D}\left(\vec{r},t\right) = \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \vec{E}\left(\vec{r},t\right) + \varepsilon_0 \chi^3\left(\vec{r}\right)\left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right)$$ 

$$\vec{B}\left(\vec{r},t\right) = \mu_0 \vec{H}\left(\vec{r},t\right)$$

We then insert these constitutive relations in the maxwell equations and arrive at the system:

$$\nabla \cdot \left( \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \vec{E}\left(\vec{r},t\right) + \varepsilon_0 \chi^3\left(\vec{r}\right)\left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right) \right) = 0$$

$$\nabla \times \vec{E}\left(\vec{r},t\right) = - \mu_0 \frac{\partial}{\partial t}\vec{H}\left(\vec{r},t\right)$$

$$\nabla \times \vec{H}\left(\vec{r},t\right) = \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \frac{\partial}{\partial t} \vec{E}\left(\vec{r},t\right) + \varepsilon_0 \chi^3\left(\vec{r}\right) \frac{\partial}{\partial t} \left( \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right)\right) + \sigma\left(\vec{r},t\right) \vec{E}\left(\vec{r},t\right)$$

$$\nabla \cdot \vec{H}\left(\vec{r},t\right) = 0$$

We have to simplify to get a solvable form:

$$\nabla \times \nabla \times \vec{E}\left(\vec{r},t\right) = - \mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \frac{\partial^2}{\partial t^2} \vec{E}\left(\vec{r},t\right) - \mu_0 \varepsilon_0 \chi^3\left(\vec{r}\right) \frac{\partial^2}{\partial t^2} \left( \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right)\right) - \mu_0 \frac{\partial \sigma\left(\vec{r},t\right)}{\partial t} \vec{E}\left(\vec{r},t\right)- \mu_0 \frac{\partial \vec{E}\left(\vec{r},t\right)}{\partial t}\sigma\left(\vec{r},t\right)$$

$$\nabla \times \nabla \times \vec{E}\left(\vec{r},t\right) = \nabla \left(\nabla \cdot \vec{E}\left(\vec{r},t\right)\right) - \nabla^2 \vec{E}\left(\vec{r},t\right)$$

$$\nabla \cdot \vec{E}\left(\vec{r},t\right) = - \frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} - \frac{1}{\varepsilon \left(\vec{r}\right)} \nabla \cdot \left(\chi^3\left(\vec{r}\right)\left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right) \right)$$

$$\nabla \left(\nabla \cdot \vec{E}\left(\vec{r},t\right)\right) - \nabla^2 \vec{E}\left(\vec{r},t\right) = - \nabla^2 \vec{E}\left(\vec{r},t\right) + \nabla \left( - \frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} - \frac{1}{\varepsilon \left(\vec{r}\right)} \nabla \cdot \left(\chi^3\left(\vec{r}\right)\left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right) \right) \right)$$


$$\frac{\partial^2}{\partial t^2} \left( \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right)\right) = \frac{\partial}{\partial t} \left( \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \frac{\partial}{\partial t} \vec{E}\left(\vec{r},t\right) + 2 \left(\sum \frac{\partial}{\partial t} \left\lvert E_i \left(\vec{r},t\right) \right\rvert \right)\vec{E}\left(\vec{r},t\right)\right)$$


$$\frac{\partial^2}{\partial t^2} \left( \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right)\right) = \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \frac{\partial^2}{\partial t^2} \vec{E}\left(\vec{r},t\right) + 4 \left(\sum \frac{\partial}{\partial t} \left\lvert E_i \left(\vec{r},t\right) \right\rvert \right) \frac{\partial}{\partial t}\vec{E}\left(\vec{r},t\right) + 2\left(\sum \frac{\partial^2}{\partial t^2} \left\lvert E_i \left(\vec{r},t\right) \right\rvert \right)\vec{E}\left(\vec{r},t\right)$$

Through the addition of the nonlinear term, we get a nasty term that makes it impossible to write the equation system in a clean matrix form as there is a sum over second-order time derivatives. We introduce an approximation to mitigate this problem. We assume that the absolute square of the electric field is a slowly varying function of time:

$$\frac{\partial^2}{\partial t^2} \left( \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right)\right) \approx \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \frac{\partial^2}{\partial t^2} \vec{E}\left(\vec{r},t\right)$$

Then we can plug everythin in one equation:

$$- \nabla^2 \vec{E}\left(\vec{r},t\right) + \nabla \left( - \frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} - \frac{1}{\varepsilon \left(\vec{r}\right)} \nabla \cdot \left(\chi^3\left(\vec{r}\right)\left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right) \right) \right) = - \mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \frac{\partial^2}{\partial t^2} \vec{E}\left(\vec{r},t\right) - \mu_0 \varepsilon_0 \chi^3\left(\vec{r}\right) \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \frac{\partial^2}{\partial t^2} \vec{E}\left(\vec{r},t\right) - \mu_0 \frac{\partial \sigma\left(\vec{r},t\right)}{\partial t} \vec{E}\left(\vec{r},t\right)- \mu_0 \frac{\partial \vec{E}\left(\vec{r},t\right)}{\partial t}\sigma\left(\vec{r},t\right)$$

As in the previous linear problem, we can write a system of equations:

$$\frac{\partial}{\partial t}\begin{bmatrix}
        \vec{u}\left(\vec{r},t\right)\\
        \vec{v}\left(\vec{r},t\right)\\
    \end{bmatrix} 
    = \begin{bmatrix}
        0 & 1\\
        -\frac{\frac{\partial}{\partial t}\sigma\left(\vec{r},t\right)}{\varepsilon_0 \varepsilon_r\left(\vec{r}\right) + \varepsilon_0 \chi^3\left(\vec{r}\right) \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2} & -\frac{\sigma\left(\vec{r},t\right)}{\varepsilon_0 \varepsilon_r\left(\vec{r}\right) + \varepsilon_0 \chi^3\left(\vec{r}\right) \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2} \\
    \end{bmatrix}
    \begin{bmatrix}
        \vec{u}\left(\vec{r},t\right)\\
        \vec{v}\left(\vec{r},t\right)\\
    \end{bmatrix} + 
    \begin{bmatrix}
        0\\
        \frac{1}{\mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) + \mu_0\varepsilon_0 \chi^3\left(\vec{r}\right) \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2}\nabla^2 \vec{u}\left(\vec{r},t\right) + \frac{1}{\mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) + \mu_0\varepsilon_0 \chi^3\left(\vec{r}\right) \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2} \nabla \left(\frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{u}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} \right) + \frac{1}{\mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) + \mu_0\varepsilon_0 \chi^3\left(\vec{r}\right) \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2} \nabla \left(\frac{1}{\varepsilon \left(\vec{r}\right)} \nabla \cdot \left(\chi^3\left(\vec{r}\right)\left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right) \right) \right)\\
    \end{bmatrix}$$

$$\vec{u}\left(\vec{r},t\right) = \vec{E}\left(\vec{r},t\right)$$

$$\vec{v}\left(\vec{r},t\right) = \frac{\partial}{\partial t} \vec{E}\left(\vec{r},t\right)$$

# Boundary Condition

In this project, we use various types of boundary condition for modeling different physical systems:

## Perfect Electric Conductor Boundary Conditions
This boundary condition models a perfect metallic mirror. In The whole field is reflected at the boundary and no field is allowed to exist inside the mirror. The perfect electric conductor is represented through the dirichlet boundary condition:

$$\vec{E}\left(\vec{r},t\right) = 0, \; \vec{r} \in \partial \Omega_{PEC}$$

## Perfect Absorbing Boundary Conditions

Perfect absorbing boundary conditions should model infinit free space propagation of the field without reflections.

A first naive idea is to use again the perfect electric conductor boundary conditions, but introduce a highly absorbing conductivity layer at the boundary.

Another choice for the perfect absorbing boundary conditions is the zeroth order approximation. This choice is represented through the following mixed boundary condition:

$$\frac{\partial}{\partial t}\vec{E}\left(\vec{r},t\right) + c \left(\vec{n} \cdot \nabla \right)\vec{E}\left(\vec{r},t\right), \; \vec{r} \in \partial \Omega_{PA0}$$

Where $c$ is the speed of light in the boundary medium and $\vec{n}$ is the outwards pointing boundary normal vector. This can be then rewriten in the following way:

$$\vec{v}\left(\vec{r},t\right) = - c \left(\vec{n} \cdot \nabla \right)\vec{u}\left(\vec{r},t\right), \; \vec{r} \in \partial \Omega_{PA0}$$

Therefore, the velocity can be directly calculated on the boundary and is not undefined anymore as before. In addition, we need the value of $\vec{u}\left(\vec{r},t\right)$ on the boundary. This we get through applying the update rule not only in the inside of the domain, but on the absorbing boundary.

# Results
In the following chapter, we present plots and gifs resulting from running different physics examples. All of the 2D examples were run on four GPUs, whereas the 3D examples were run on eight GPUs:

## Homogenous 2D Problem
For the first example, resulting from running:
```bash
julia --project examples/homogeneous_dirichlet_Ez_2D.jl
julia --project examples/homogeneous_dirichlet_Ez_2D_viz.jl
``` 
We simulated in 2D the homogenous wave equation in vacuum with a initial condition of a modulated cosine modulated gaussian pulse with only perfect reflecting boundary conditions. 

![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/homogeneous_dirichlet_Ez.gif)

## Homogenous 3D Problem
```bash
julia --project examples/homogeneous_dirichlet_vecE_3D.jl
julia --project examples/vecE_3D_viz.jl
``` 
Where one has to change the settings in vecE_3D_viz.jl to match the previous script.

We simulate the derived system of equation in the linear lossy problem section with only perfect reflecting boundary conditions. 

With the following normalized static material paramters:
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/homogeneous_dirichlet_vecE_3D_epsilon.png)
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/homogeneous_dirichlet_vecE_3D_sigma.png)

And a y-polarized field with following initial conditions:
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/homogeneous_dirichlet_vecE_3D_inituy.png)
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/homogeneous_dirichlet_vecE_3D_sliceplty.gif)

Then we present the propagation of the x- and y- component of the field since coupling to the z-component is vanishing with the choosen material parameters:
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/homogeneous_dirichlet_vecE_3D_slicepltx.gif)
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/homogeneous_dirichlet_vecE_3D_initvy.png)



## Absorbing 3D Problem
```bash
julia --project examples/absorbing_boundary_vecE_3D.jl
julia --project examples/vecE_3D_viz.jl
``` 
Where one has to change the settings in vecE_3D_viz.jl to match the previous script.

We simulate the derived system of equation in the linear lossy problem section. In contrast to the previous section, on the left yz plane a first order perfect absorbing boundary conditions is applied. 

With the following normalized static material paramters:
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/absorbing_boundary_vecE_3D_epsilon.png)
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/absorbing_boundary_vecE_3D_sigma.png)

And a y-polarized field with following initial conditions:
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/absorbing_boundary_vecE_3D_inituy.png)
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/absorbing_boundary_vecE_3D_sliceplty.gif)

Then we present the propagation of the x- and y- component of the field since coupling to the z-component is vanishing with the choosen material parameters:
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/absorbing_boundary_vecE_3D_slicepltx.gif)
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/absorbing_boundary_vecE_3D_initvy.png)

## Nonlinear 3D Problem
```bash
julia --project examples/nonlinear_vecE_3D.jl
julia --project examples/vecE_3D_viz.jl
``` 
Where one has to change the settings in vecE_3D_viz.jl to match the previous script.

We simulate the derived system of equation in the $\chi^3$ nonlinear lossy problem section with only perfect reflecting boundary conditions. 

With the following normalized static material paramters:
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/nonlinear_vecE_3D_epsilon.png)
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/nonlinear_vecE_3D_sigma.png)
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/nonlinear_vecE_3D_chi3.png)


And a y-polarized field with following initial conditions:
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/nonlinear_vecE_3D_inituy.png)
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/nonlinear_vecE_3D_sliceplty.gif)

Then we present the propagation of the x- and y- component of the field since coupling to the z-component is vanishing with the choosen material parameters:
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/nonlinear_vecE_3D_slicepltx.gif)
![](https://github.com/lheimabch/MaederHeimbach/blob/main/MaxwellWave/docs/nonlinear_vecE_3D_initvy.png)

# Conclusion and Possible Extensions

# Running the Project Examples
