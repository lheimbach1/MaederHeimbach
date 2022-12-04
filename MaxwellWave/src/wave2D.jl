@doc raw"""
    @parallel_indices (ix,iy) function update_u(u,v,dt)

Computes the update rule for the scalar field through the equation:
```math
u^{j+1} = u^{j} + dt*v^{j+\frac{1}{2}}
```
"""
@parallel_indices (ix,iy) function update_u!(u,v,dt)
    nx,ny=size(u) 
    if(ix>=2 && iy>=2 && ix<=nx-1 && iy<=ny-1)
        u[ix,iy] = u[ix,iy] + dt*v[ix,iy] 
    end
    return nothing
end
@doc raw"""
    @parallel_indices (ix,iy) function update_v!(u,v,dt,_dx,_dy,c2)

Computes the update rule for the velocity of the scalar field through the equation:
```math
v^{j+\frac{1}{2}} = v^{j-\frac{1}{2}} + dt*c^2 \nabla^2 u^{j}
```
"""
@parallel_indices (ix,iy) function update_v!(u,v,dt,_dx,_dy,c2)
    nx,ny=size(u) 
    if(ix>=2 && iy>=2 && ix<=nx-1 && iy<=ny-1)
        v[ix,iy] = v[ix,iy] + dt*c2*( (u[ix+1,iy] -2*u[ix,iy] + u[ix-1,iy])*_dx + (u[ix,iy+1] - 2*u[ix,iy] + u[ix,iy-1])*_dy) 
    end
    return nothing
end