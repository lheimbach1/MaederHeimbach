# Maxwell Wave
The goal of this project is to solve the electromagnetic wave equation in 3D.

# Derivation
We start with the general Maxwell's equation:<br /> 

$$\nabla \cdot \vec{D}\left(\vec{r},t\right) = p_0\left(\vec{r},t\right)$$

$$\nabla \times \vec{E}\left(\vec{r},t\right) = - \frac{\partial}{\partial t}\vec{B}\left(\vec{r},t\right)$$

$$\nabla \times \vec{H}\left(\vec{r},t\right) = \frac{\partial}{\partial t} \vec{D}\left(\vec{r},t\right) + \vec{j}_0\left(\vec{r},t\right)$$

$$\nabla \cdot \vec{B}\left(\vec{r},t\right) = 0$$

We have to further supplement this system of equation with constitutive relations.<br /> We choose for the first problem setting a nonmagnetic dielectric with no free charges and no nonlinearities:<br />

$$p_0\left(\vec{r},t\right) = 0$$

$$\vec{j}_0\left(\vec{r},t\right) = \sigma\left(\vec{r},t\right) \vec{E}\left(\vec{r},t\right)$$

$$\vec{D}\left(\vec{r},t\right) = \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \vec{E}\left(\vec{r},t\right)$$

$$\vec{B}\left(\vec{r},t\right) = \mu_0 \vec{H}\left(\vec{r},t\right)$$

We then insert these constitutive relations in the maxwell equations and arrive at the system:<br /> 

$$\nabla \cdot \left( \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \vec{E}\left(\vec{r},t\right) \right) = 0$$

$$\nabla \times \vec{E}\left(\vec{r},t\right) = - \mu_0 \frac{\partial}{\partial t}\vec{H}\left(\vec{r},t\right)$$

$$\nabla \times \vec{H}\left(\vec{r},t\right) = \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \frac{\partial}{\partial t} \vec{E}\left(\vec{r},t\right) + \sigma\left(\vec{r},t\right) \vec{E}\left(\vec{r},t\right)$$

$$\nabla \cdot \vec{H}\left(\vec{r},t\right) = 0$$

To further simplify the system, we apply multiple transformations and vector caculus identities:<br /> 

$$\nabla \times \nabla \times \vec{E}\left(\vec{r},t\right) = - \mu_0 \frac{\partial}{\partial t} \nabla \times \vec{H}\left(\vec{r},t\right)$$

$$\nabla \times \nabla \times \vec{E}\left(\vec{r},t\right) = - \mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \frac{\partial^2}{\partial t^2}\vec{E}\left(\vec{r},t\right) -\mu_0 \left(\frac{\partial}{\partial t}\sigma\left(\vec{r},t\right)\right)\vec{E}\left(\vec{r},t\right) -\mu_0\sigma\left(\vec{r},t\right) \left(\frac{\partial}{\partial t}\vec{E}\left(\vec{r},t\right)\right)$$

$$\nabla \times \nabla \times \vec{E}\left(\vec{r},t\right) = \nabla \left(\nabla \cdot \vec{E}\left(\vec{r},t\right)\right) - \nabla^2 \vec{E}\left(\vec{r},t\right)$$

$$0 = \nabla \cdot \left(\varepsilon_r\left(\vec{r}\right) \vec{E}\left(\vec{r},t\right) \right) = \varepsilon_r\left(\vec{r}\right) \nabla \cdot \vec{E}\left(\vec{r},t\right) + \nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)$$

$$\nabla \cdot \vec{E}\left(\vec{r},t\right) = \frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)}$$

$$\nabla \times \nabla \times \vec{E}\left(\vec{r},t\right) = \nabla \left(\frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} \right) - \nabla^2 \vec{E}\left(\vec{r},t\right)$$

With all of these simplifications, we arrive at the final equation for the electric field, which we want to solve: <br /> 

$$\nabla^2 \vec{E}\left(\vec{r},t\right) - \nabla \left(\frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} \right) - \mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \frac{\partial^2}{\partial t^2}\vec{E}\left(\vec{r},t\right) -\mu_0 \left(\frac{\partial}{\partial t}\sigma\left(\vec{r},t\right)\right)\vec{E}\left(\vec{r},t\right) -\mu_0\sigma\left(\vec{r},t\right) \left(\frac{\partial}{\partial t}\vec{E}\left(\vec{r},t\right)\right) = 0$$

which can be equivalently formulated as: <br /> 

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
        \frac{1}{\mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right)}\nabla^2 \vec{u}\left(\vec{r},t\right) - \frac{1}{\mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right)} \nabla \left(\frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{u}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} \right)\\
    \end{bmatrix}$$

$$\vec{u}\left(\vec{r},t\right) = \vec{E}\left(\vec{r},t\right)$$

$$\vec{v}\left(\vec{r},t\right) = \frac{\partial}{\partial t} \vec{E}\left(\vec{r},t\right)$$