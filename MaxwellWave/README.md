# Maxwell Wave
The goal of this project is to solve the electromagnetic wave equation in 3D.

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

As in the previous nonlinear problem, we can write a system of equations:

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

Another choice for the boundary of the domain are perfect absorbing boundary conditions. These models infinit free space propagation of the field without reflections. The following mixed boundary condition is a possible realisation:

$$\vec{E}\left(\vec{r},t\right) = todo, \; \vec{r} \in \partial \Omega_{PA}$$

# Results

# Conclusion and Possible Extensions
