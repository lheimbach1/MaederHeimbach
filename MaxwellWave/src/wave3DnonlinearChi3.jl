@doc raw"""
    @parallel_indices (ix, iy, iz) function update_vecv_chi3!(ux, uy, uz, vx, vy, vz, dt, _dx, _dy, _dz, _dx_2, _dy_2, _dz_2, alpha, delta, omega)

Computes the update rule for the velocity of the vector field through the equation:
```math
\vec{v}^{j+\frac{1}{2}} = \vec{v}^{j-\frac{1}{2}} + dt*\alpha\left(\vec{r}\right) \nabla \left(\delta \left(\vec{r}\right) \nabla \cdot \left( \omega \left(\vec{r}\right) \left\lvert \vec{u}^{j} \right\rvert^2 \vec{u}^{j} \right) \right)
```

and discretized in the following way:

```math
v_x^{j+\frac{1}{2}}\left[ix,iy,iz\right] = v_x^{j-\frac{1}{2}}\left[ix,iy,iz\right] + dt*\alpha\left[ix,iy,iz\right] \left( 
\delta \left[ix+1,iy,iz\right] \left( \left(\omega \left[ix+1,iy,iz\right] \left\lvert \vec{u}^{j}\left[ix+1,iy,iz\right] \right\rvert^2 u_x^j \left[ix+1,iy,iz\right]  
 - \omega \left[ix,iy,iz\right] \left\lvert \vec{u}^{j}\left[ix,iy,iz\right] \right\rvert^2 u_x^j\left[ix,iy,iz\right]\right)/dx + \left(\omega \left[ix+1,iy+1,iz\right] 
 \left\lvert \vec{u}^{j}\left[ix+1,iy+1,iz\right] \right\rvert^2 u_y^j\left[ix+1,iy+1,iz\right]  
 - \omega \left[ix+1,iy-1,iz\right] \left\lvert \vec{u}^{j}\left[ix+1,iy-1,iz\right] \right\rvert^2 u_y^j\left[ix+1,iy-1,iz\right]\right)/dy/2+ \left(\omega \left[ix+1,iy,iz+1\right] \left\lvert \vec{u}^{j}\left[ix+1,iy,iz+1\right] \right\rvert^2 u_z^j\left[ix+1,iy,iz+1\right]- \omega \left[ix+1,iy,iz-1\right] \left\lvert \vec{u}^{j}\left[ix+1,iy,iz-1\right] \right\rvert^2 u_z^j\left[ix+1,iy,iz-1\right]\right)/dz/2\right) + \delta \left[ix-1,iy,iz\right] \left( \left(\omega \left[ix,iy,iz\right] \left\lvert \vec{u}^{j}\left[ix,iy,iz\right] \right\rvert^2 u_x^j\left[ix,iy,iz\right] 
 - \omega \left[ix-1,iy,iz\right] \left\lvert \vec{u}^{j}\left[ix-1,iy,iz\right] \right\rvert^2 u_x^j\left[ix-1,iy,iz\right]\right)/dx 
+ \left(\omega \left[ix-1,iy+1,iz\right] \left\lvert \vec{u}^{j}\left[ix-1,iy+1,iz\right] \right\rvert^2 u_y^j\left[ix-1,iy+1,iz\right]  
 - \omega \left[ix-1,iy-1,iz\right] \left\lvert \vec{u}^{j}\left[ix-1,iy-1,iz\right] \right\rvert^2 u_y^j\left[ix-1,iy-1,iz\right]\right)/dy/2 + \left(\omega \left[ix-1,iy,iz+1\right] \left\lvert \vec{u}^{j}\left[ix-1,iy,iz+1\right] \right\rvert^2 u_z^j\left[ix-1,iy,iz+1\right]- \omega \left[ix-1,iy,iz-1\right] \left\lvert \vec{u}^{j}\left[ix-1,iy,iz-1\right] \right\rvert^2 u_z^j\left[ix-1,iy,iz-1\right]\right)/dz/2\right) \right)/dx/2
```

```math
v_y^{j+\frac{1}{2}}\left[ix,iy,iz\right] = v_y^{j-\frac{1}{2}}\left[ix,iy,iz\right] + dt*\alpha\left[ix,iy,iz\right] \left( 
\delta \left[ix,iy+1,iz\right] \left( \left(\omega \left[ix+1,iy+1,iz\right] \left\lvert \vec{u}^{j}\left[ix+1,iy+1,iz\right] \right\rvert^2 u_x^j \left[ix+1,iy+1,iz\right]  
 - \omega \left[ix-1,iy+1,iz\right] \left\lvert \vec{u}^{j}\left[ix-1,iy+1,iz\right] \right\rvert^2 u_x^j\left[ix-1,iy+1,iz\right]\right)/dx/2 + \left(\omega \left[ix,iy+1,iz\right] 
 \left\lvert \vec{u}^{j}\left[ix,iy+1,iz\right] \right\rvert^2 u_y^j\left[ix,iy+1,iz\right]  
 - \omega \left[ix,iy,iz\right] \left\lvert \vec{u}^{j}\left[ix,iy,iz\right] \right\rvert^2 u_y^j\left[ix,iy,iz\right]\right)/dy+ 
 \left(\omega \left[ix,iy+1,iz+1\right] \left\lvert \vec{u}^{j}\left[ix,iy+1,iz+1\right] \right\rvert^2 u_z^j\left[ix,iy+1,iz+1\right]- \omega \left[ix,iy+1,iz-1\right] \left\lvert \vec{u}^{j}\left[ix,iy+1,iz-1\right] \right\rvert^2 u_z^j\left[ix,iy+1,iz-1\right]\right)/dz/2\right) + \delta \left[ix,iy-1,iz\right] \left( \left(\omega \left[ix+1,iy-1,iz\right] \left\lvert \vec{u}^{j}\left[ix+1,iy-1,iz\right] \right\rvert^2 u_x^j\left[ix+1,iy-1,iz\right] 
 - \omega \left[ix-1,iy-1,iz\right] \left\lvert \vec{u}^{j}\left[ix-1,iy-1,iz\right] \right\rvert^2 u_x^j\left[ix-1,iy-1,iz\right]\right)/dx/2 
+ \left(\omega \left[ix,iy,iz\right] \left\lvert \vec{u}^{j}\left[ix,iy,iz\right] \right\rvert^2 u_y^j\left[ix,iy,iz\right]  
 - \omega \left[ix,iy-1,iz\right] \left\lvert \vec{u}^{j}\left[ix,iy-1,iz\right] \right\rvert^2 u_y^j\left[ix,iy-1,iz\right]\right)/dy + \left(\omega \left[ix,iy-1,iz+1\right] \left\lvert \vec{u}^{j}\left[ix,iy-1,iz+1\right] \right\rvert^2 u_z^j\left[ix,iy-1,iz+1\right]- \omega \left[ix,iy-1,iz-1\right] \left\lvert \vec{u}^{j}\left[ix,iy-1,iz-1\right] \right\rvert^2 u_z^j\left[ix,iy-1,iz-1\right]\right)/dz/2\right) \right)/dy/2
```

```math
v_z^{j+\frac{1}{2}}\left[ix,iy,iz\right] = v_z^{j-\frac{1}{2}}\left[ix,iy,iz\right] + dt*\alpha\left[ix,iy,iz\right] \left( 
\delta \left[ix,iy,iz+1\right] \left( \left(\omega \left[ix+1,iy,iz+1\right] \left\lvert \vec{u}^{j}\left[ix+1,iy,iz+1\right] \right\rvert^2 u_x^j \left[ix+1,iy,iz+1\right]  
 - \omega \left[ix-1,iy,iz+1\right] \left\lvert \vec{u}^{j}\left[ix-1,iy,iz+1\right] \right\rvert^2 u_x^j\left[ix-1,iy,iz+1\right]\right)/dx/2 + \left(\omega \left[ix,iy+1,iz+1\right] 
 \left\lvert \vec{u}^{j}\left[ix,iy+1,iz+1\right] \right\rvert^2 u_y^j\left[ix,iy+1,iz+1\right]  
 - \omega \left[ix,iy-1,iz+1\right] \left\lvert \vec{u}^{j}\left[ix,iy-1,iz+1\right] \right\rvert^2 u_y^j\left[ix,iy-1,iz+1\right]\right)/dy/2+ 
 \left(\omega \left[ix,iy,iz+1\right] \left\lvert \vec{u}^{j}\left[ix,iy,iz+1\right] \right\rvert^2 u_z^j\left[ix,iy,iz+1\right]- 
 \omega \left[ix,iy,iz\right] \left\lvert \vec{u}^{j}\left[ix,iy,iz\right] \right\rvert^2 u_z^j\left[ix,iy,iz\right]\right)/dz\right) +
 \delta \left[ix,iy,iz-1\right] \left( \left(\omega \left[ix+1,iy,iz-1\right] \left\lvert \vec{u}^{j}\left[ix+1,iy,iz-1\right] \right\rvert^2 u_x^j\left[ix+1,iy,iz-1\right] 
 - \omega \left[ix-1,iy,iz-1\right] \left\lvert \vec{u}^{j}\left[ix-1,iy,iz-1\right] \right\rvert^2 u_x^j\left[ix-1,iy,iz-1\right]\right)/dx/2 
+ \left(\omega \left[ix,iy+1,iz-1\right] \left\lvert \vec{u}^{j}\left[ix,iy+1,iz-1\right] \right\rvert^2 u_y^j\left[ix,iy+1,iz-1\right]  
 - \omega \left[ix,iy-1,iz-1\right] \left\lvert \vec{u}^{j}\left[ix,iy-1,iz-1\right] \right\rvert^2 u_y^j\left[ix,iy-1,iz-1\right]\right)/dy/2 + \left(\omega \left[ix,iy,iz\right] \left\lvert \vec{u}^{j}\left[ix,iy,iz\right] \right\rvert^2 u_z^j\left[ix,iy,iz\right]- \omega 
 \left[ix,iy,iz-1\right] \left\lvert \vec{u}^{j}\left[ix,iy,iz-1\right] \right\rvert^2 u_z^j\left[ix,iy,iz-1\right]\right)/dz\right) \right)/dz/2
```

Where we have applied a central difference formula for the outer derivative. For the inner derivates, we use the central or forward/backward formula depending on if the derivative is in the same direction as the outer derivative. We use the forward/backward formula to arrive at a local nearest neighbor stencil.

"""
@parallel_indices (ix, iy, iz) function update_vecv_chi3!(ux, uy, uz, vx, vy, vz, dt, _dx, _dy, _dz, _dx_2, _dy_2, _dz_2, alpha, delta, omega)
    nx, ny, nz = size(ux)
    if (ix >= 2 && iy >= 2 && iz >= 2 && ix <= nx - 1 && iy <= ny - 1 && iz <= nz - 1)

        # load the nearest neighbor values into local memory
        # using mutable static arrays as a workaround the problem of memory allocation in the gpu kernel
        ux_local = MArray{Tuple{3,3,3},Float64,3,27}(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        uy_local = MArray{Tuple{3,3,3},Float64,3,27}(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        uz_local = MArray{Tuple{3,3,3},Float64,3,27}(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)

        # alocate memory for absolute value
        abs2_u_local = MArray{Tuple{3,3,3},Float64,3,27}(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)

        # compute nearest neighbor absolute values 
        for ix_locla = 1:3
            for iy_locla = 1:3
                for iz_locla = 1:3
                    ux_local[ix_locla,iy_locla,iz_locla] = ux[ix-2+ix_locla,iy-2+iy_locla,iz-2+iz_locla]
                    uy_local[ix_locla,iy_locla,iz_locla] = uy[ix-2+ix_locla,iy-2+iy_locla,iz-2+iz_locla]
                    uz_local[ix_locla,iy_locla,iz_locla] = uz[ix-2+ix_locla,iy-2+iy_locla,iz-2+iz_locla]
                    abs2_u_local[ix_locla,iy_locla,iz_locla] = ux_local[ix_locla,iy_locla,iz_locla]*ux_local[ix_locla,iy_locla,iz_locla] + 
                                                               uy_local[ix_locla,iy_locla,iz_locla]*uy_local[ix_locla,iy_locla,iz_locla] +
                                                               uz_local[ix_locla,iy_locla,iz_locla]*uz_local[ix_locla,iy_locla,iz_locla]
                end 
            end 
        end    

        # implement the update rule
        vx[ix, iy, iz] = vx[ix, iy, iz] + dt * alpha[ix, iy, iz] * (delta[ix+1,iy,iz]*
                                                                    ( (omega[ix+1,iy,iz]*abs2_u_local[3,2,2]*ux_local[3,2,2] -
                                                                       omega[ix,iy,iz]*abs2_u_local[2,2,2]*ux_local[2,2,2])*_dx +
                                                                      (omega[ix+1,iy+1,iz]*abs2_u_local[3,3,2]*uy_local[3,3,2]-
                                                                       omega[ix+1,iy-1,iz]*abs2_u_local[3,1,2]*uy_local[3,1,2])*_dy_2 +
                                                                      (omega[ix+1,iy,iz+1]*abs2_u_local[3,2,3]*uz_local[3,2,3]-
                                                                       omega[ix+1,iy,iz-1]*abs2_u_local[3,2,1]*uz_local[3,2,1])*_dz_2 ) - 
                                                                    delta[ix-1,iy,iz]*
                                                                    ( (omega[ix,iy,iz]*abs2_u_local[2,2,2]*ux_local[2,2,2] -
                                                                       omega[ix-1,iy,iz]*abs2_u_local[1,2,2]*ux_local[1,2,2])*_dx +
                                                                      (omega[ix-1,iy+1,iz]*abs2_u_local[1,3,2]*uy_local[1,3,2]-
                                                                       omega[ix-1,iy-1,iz]*abs2_u_local[1,1,2]*uy_local[1,1,2])*_dy_2 +
                                                                      (omega[ix-1,iy,iz+1]*abs2_u_local[1,2,3]*uz_local[1,2,3]-
                                                                       omega[ix-1,iy,iz-1]*abs2_u_local[1,2,1]*uz_local[1,2,1])*_dz_2 ) )*_dx_2

        vy[ix, iy, iz] = vy[ix, iy, iz] + dt * alpha[ix, iy, iz] * (delta[ix,iy+1,iz]*
                                                                    ( (omega[ix+1,iy+1,iz]*abs2_u_local[3,3,2]*ux_local[3,3,2] -
                                                                       omega[ix-1,iy+1,iz]*abs2_u_local[1,3,2]*ux_local[1,3,2])*_dx_2 +
                                                                      (omega[ix,iy+1,iz]*abs2_u_local[2,3,2]*uy_local[2,3,2]-
                                                                       omega[ix,iy,iz]*abs2_u_local[2,2,2]*uy_local[2,2,2])*_dy +
                                                                      (omega[ix,iy+1,iz+1]*abs2_u_local[2,3,3]*uz_local[2,3,3]-
                                                                       omega[ix,iy+1,iz-1]*abs2_u_local[2,3,1]*uz_local[2,3,1])*_dz_2 ) - 
                                                                    delta[ix,iy-1,iz]*
                                                                    ( (omega[ix+1,iy-1,iz]*abs2_u_local[3,1,2]*ux_local[3,1,2] -
                                                                       omega[ix-1,iy-1,iz]*abs2_u_local[1,1,2]*ux_local[1,1,2])*_dx_2 +
                                                                      (omega[ix,iy,iz]*abs2_u_local[2,2,2]*uy_local[2,2,2]-
                                                                       omega[ix,iy-1,iz]*abs2_u_local[2,1,2]*uy_local[2,1,2])*_dy +
                                                                      (omega[ix,iy-1,iz+1]*abs2_u_local[2,1,3]*uz_local[2,1,3]-
                                                                       omega[ix,iy-1,iz-1]*abs2_u_local[2,1,1]*uz_local[2,1,1])*_dz_2 ) )*_dy_2


        vz[ix, iy, iz] = vz[ix, iy, iz] + dt * alpha[ix, iy, iz] * (delta[ix,iy,iz+1]*
                                                                    ( (omega[ix+1,iy,iz+1]*abs2_u_local[3,2,3]*ux_local[3,2,3] -
                                                                       omega[ix-1,iy,iz+1]*abs2_u_local[1,2,3]*ux_local[1,2,3])*_dx_2 +
                                                                      (omega[ix,iy+1,iz+1]*abs2_u_local[2,3,3]*uy_local[2,3,3]-
                                                                       omega[ix,iy-1,iz+1]*abs2_u_local[2,1,3]*uy_local[2,1,3])*_dy_2 +
                                                                      (omega[ix,iy,iz+1]*abs2_u_local[2,2,3]*uz_local[2,2,3]-
                                                                       omega[ix,iy,iz]*abs2_u_local[2,2,2]*uz_local[2,2,2])*_dz ) - 
                                                                    delta[ix,iy,iz-1]*
                                                                    ( (omega[ix+1,iy,iz-1]*abs2_u_local[3,2,1]*ux_local[3,2,1] -
                                                                       omega[ix-1,iy,iz-1]*abs2_u_local[1,2,1]*ux_local[1,2,1])*_dx_2 +
                                                                      (omega[ix,iy+1,iz-1]*abs2_u_local[2,3,1]*uy_local[2,3,1]-
                                                                       omega[ix,iy-1,iz-1]*abs2_u_local[2,1,1]*uy_local[2,1,1])*_dy_2 +
                                                                      (omega[ix,iy,iz]*abs2_u_local[2,2,2]*uz_local[2,2,2]-
                                                                       omega[ix,iy,iz-1]*abs2_u_local[2,2,1]*uz_local[2,2,1])*_dz ) )*_dz_2


    end
    return nothing
end