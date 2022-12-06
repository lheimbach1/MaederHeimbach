
@doc raw"""
homogeneous_dirichlet_Ez_2D(; do_visu=false)

Solves the two wave equation problem:

```math
\nabla^2 \vec{E} - \frac{1}{c^2} \frac{\partial^2 }{\partial t^2} \vec{E} = 0, \; \vec{E} \in \Omega
```
```math
\vec{E} = 0, \; \vec{E} \in \partial \Omega
```

on a regular grid with ``c = 1`` and ``L_x = L_y = 10``

Formulated as: 
```math
\frac{\partial }{\partial t} u = v
```
```math
\frac{\partial }{\partial t} v = c^2 \nabla^2 u
```
and discretized in the leap frog fashion to ensure approximatly energy conservation:
```math
v^{j+\frac{1}{2}} = v^{j-\frac{1}{2}} + dt*c^2 \nabla^2 u^{j}
```
```math
u^{j+1} = u^{j} + dt*v^{j+\frac{1}{2}}
```

"""
function homogeneous_dirichlet_vecE_3D(; do_visu=false)

end


homogeneous_dirichlet_vecE_3D(do_visu=true)