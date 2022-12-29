# instantiate for running on daint
if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.instantiate()
end


using Documenter

using Plots, ParallelStencil, ParallelStencil.FiniteDifferences3D, Printf
using ImplicitGlobalGrid
using StaticArrays
import MPI
ENV["GKSwstype"] = "nul"

# decide if one uses the gpu or cpu
const USE_GPU = false
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3)
else
    @init_parallel_stencil(Threads, Float64, 3)
end


include("../src/wave3D.jl")
include("../src/wave3DnonlinearChi3.jl")
include("../src/auxiliary.jl")

@doc raw"""
    nonlinear_vecE_3D(; do_visu=false)
Solves the 3D wave equation problem:
```math
- \nabla^2 \vec{E}\left(\vec{r},t\right) + \nabla \left( - \frac{\nabla \varepsilon_r\left(\vec{r}\right) \cdot \vec{E}\left(\vec{r},t\right)}{\varepsilon_r\left(\vec{r}\right)} - \frac{1}{\varepsilon \left(\vec{r}\right)} \nabla \cdot \left(\chi^3\left(\vec{r}\right)\left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \vec{E}\left(\vec{r},t\right) \right) \right) = - \mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) \frac{\partial^2}{\partial t^2} \vec{E}\left(\vec{r},t\right) - \mu_0 \varepsilon_0 \chi^3\left(\vec{r}\right) \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2 \frac{\partial^2}{\partial t^2} \vec{E}\left(\vec{r},t\right) - \mu_0 \frac{\partial \sigma\left(\vec{r},t\right)}{\partial t} \vec{E}\left(\vec{r},t\right)- \mu_0 \frac{\partial \vec{E}\left(\vec{r},t\right)}{\partial t}\sigma\left(\vec{r},t\right)
```
```math
\vec{E} = 0, \; \vec{E} \in \partial \Omega \times \left[0,T\right]
```
```math
\vec{E} = 2*\vec{n}_y e^{-(y - L_y/2)^2 -(z - L_z/2)^2}  sech\left(0.5 \left(k_x x - \omega t \right)\right), \; E_z \in  \Omega \times \left[0\right]
```
on a regular grid with ``\varepsilon_0 = 1``, ``\varepsilon_r\left(\vec{r}\right) = e^{-(y-L_y/2)^2-(z-L_z/2)^2}``, ``\sigma\left(t,\vec{r}\right) = e^{-(y-L_y/2)^2-(z-L_z/2)^2}``   and ``L_x = L_y = L_z = 10``
Formulated as: 
```math
\frac{\partial}{\partial t}\begin{bmatrix}
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
    \end{bmatrix}
```
and discretized in the leap frog fashion to ensure approximatly energy conservation:
```math
\vec{v}^{j+\frac{1}{2}} = \vec{v}^{j-\frac{1}{2}} + \vec{f}\left(\vec{u}^{j},\vec{v}^{j-\frac{1}{2}}\right)  
```
```math
\vec{u}^{j+1} = \vec{u}^{j} + dt*\vec{v}^{j+\frac{1}{2}}
```
"""
function nonlinear_vecE_3D(; do_visu=false)

    # where to save and example name
    example_name = "nonlinear_vecE_3D"

    # include the physics and numerics settings
    include(example_name * "_settings.jl")

    # initailize global grid and such MPI
    me, dims, _, _, comm_cart = init_global_grid(nx, ny, nz)
    neighbors_x = MPI.Cart_shift(comm_cart, 0, 1)
    neighbors_y = MPI.Cart_shift(comm_cart, 1, 1)
    neighbors_z = MPI.Cart_shift(comm_cart, 2, 1)

    # optimization for communication
    b_width = (8, 8, 4)

    # derived numerics
    dx = lx / nx_g()
    dy = ly / ny_g()
    dz = lz / nz_g()
    _dx2 = (1 / dx)^2
    _dy2 = (1 / dy)^2
    _dz2 = (1 / dz)^2
    _dx = 1 / dx
    _dy = 1 / dy
    _dz = 1 / dz
    _dx_2 = 1 / dx / 2
    _dy_2 = 1 / dy / 2
    _dz_2 = 1 / dz / 2
    _mu0 = 1 / mu0
    _epsilon0 = 1 / epsilon0
    dt = min(dx, dy, dy) / sqrt(c2) / 8

    # array initialisation
    ux = @zeros(nx, ny, nz)
    uy = @zeros(nx, ny, nz)
    uz = @zeros(nx, ny, nz)
    # u_pulse_shape = [cos(k0 * x_g(ix, dx, ux[:, :, :])) * exp(-(x_g(ix, dx, ux[:, :, :]) - lx / 2)^2 / (2 * sigma2[1]) - (y_g(iy, dy, ux[:, :, :]) - ly / 2)^2 / (2 * sigma2[2]) - (z_g(iz, dz, ux[:, :, :]) - lz / 2)^2 / (2 * sigma2[3])) for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)]
    u_pulse_shape = [2 * sech(0.5 * (k0 * (x_g(ix, dx, ux[:, :, :]) - lx / 2))) * exp(-(y_g(iy, dy, ux[:, :, :]) - ly / 2)^2 / (2 * sigma2[2]) - (z_g(iz, dz, ux[:, :, :]) - lz / 2)^2 / (2 * sigma2[3])) for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)]
    ux .= Data.Array(p0[1] .* u_pulse_shape)
    uy .= Data.Array(p0[2] .* u_pulse_shape)
    uz .= Data.Array(p0[3] .* u_pulse_shape)


    # set boundaries to zeros if the neighbouring rank is null
    if neighbors_x[1] == MPI.MPI_PROC_NULL
        ux[1, :, :] .= @zeros(ny, nz)
        uy[1, :, :] .= @zeros(ny, nz) 
        uz[1, :, :] .= @zeros(ny, nz)    
    end
    if neighbors_x[2] == MPI.MPI_PROC_NULL
        ux[end, :, :] .= @zeros(ny, nz)
        uy[end, :, :] .= @zeros(ny, nz)
        uz[end, :, :] .= @zeros(ny, nz)
    end
    if neighbors_y[1] == MPI.MPI_PROC_NULL
        ux[:, 1, :] .= @zeros(nx, nz)
        uy[:, 1, :] .= @zeros(nx, nz)
        uz[:, 1, :] .= @zeros(nx, nz)
    end
    if neighbors_y[2] == MPI.MPI_PROC_NULL
        ux[:, end, :] .= @zeros(nx, nz)
        uy[:, end, :] .= @zeros(nx, nz)
        uz[:, end, :] .= @zeros(nx, nz)
    end
    if neighbors_z[1] == MPI.MPI_PROC_NULL
        ux[:, :, 1] .= @zeros(nx, ny)
        uy[:, :, 1] .= @zeros(nx, ny)
        uz[:, :, 1] .= @zeros(nx, ny)
    end
    if neighbors_z[2] == MPI.MPI_PROC_NULL
        ux[:, :, end] .= @zeros(nx, ny)
        uy[:, :, end] .= @zeros(nx, ny)
        uz[:, :, end] .= @zeros(nx, ny)
    end

    # v_pulse_shape = [w0 * sin(k0 * x_g(ix, dx, ux[:, :, :])) * exp(-(x_g(ix, dx, ux[:, :, :]) - lx / 2)^2 / (2 * sigma2[1]) - (y_g(iy, dy, ux[:, :, :]) - ly / 2)^2 / (2 * sigma2[2]) - (z_g(iz, dz, ux[:, :, :]) - lz / 2)^2 / (2 * sigma2[3])) for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)]
    v_pulse_shape = [(2 * sech(0.5 * (k0 * (x_g(ix, dx, ux[:, :, :]) - lx / 2) - w0 * dt / 2)) - 2 * sech(0.5 * (k0 * (x_g(ix, dx, ux[:, :, :]) - lx / 2) + w0 * dt / 2))) / dt * exp(-(y_g(iy, dy, ux[:, :, :]) - ly / 2)^2 / (2 * sigma2[2]) - (z_g(iz, dz, ux[:, :, :]) - lz / 2)^2 / (2 * sigma2[3])) for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)]

    vx = Data.Array(p0[1] .* v_pulse_shape)
    vy = Data.Array(p0[2] .* v_pulse_shape)
    vz = Data.Array(p0[3] .* v_pulse_shape)

    # static coefficient function 

    # physics parameters
    # permitivity
    epsilon = [(1 + exp(-(y_g(iy, dy, ux[:, :, :]) - ly / 2)^2 - (z_g(iz, dz, ux[:, :, :]) - lz / 2)^2)) for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)]
    # conductivity
    sigma = [(1 - exp(-0.1*(y_g(iy, dy, ux[:, :, :]) - ly / 2)^2 - 0.1*(z_g(iz, dz, ux[:, :, :]) - lz / 2)^2)) for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)]
    # chi3 nonlinearity
    chi3 = Data.Array( [0.01*exp(-(x_g(ix, dx, ux[:, :, :]) - lx / 2)^2 -(y_g(iy, dy, ux[:, :, :]) - ly / 2)^2 - (z_g(iz, dz, ux[:, :, :]) - lz / 2)^2) for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3) ] )

    #         \frac{1}{\mu_0 \varepsilon_0 \varepsilon_r\left(\vec{r}\right) + \mu_0\varepsilon_0 \chi^3\left(\vec{r}\right) \left\lvert \vec{E}\left(\vec{r},t\right)\right\rvert^2}
    alpha = @zeros(nx, ny, nz)

    # -\frac{\frac{\partial}{\partial t}\sigma\left(\vec{r},t\right)}{\varepsilon_0 \varepsilon_r\left(\vec{r}\right)}
    beta = @zeros(nx, ny, nz)

    # -\frac{\sigma\left(\vec{r},t\right)}{\varepsilon_0 \varepsilon_r\left(\vec{r}\right)}
    gamma = Data.Array([-sigma[ix,iy,iz] * mu0  for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)])

    # \left(\frac{1}{\varepsilon \left(\vec{r}\right)}
    delta = Data.Array([1/epsilon[ix,iy,iz]  for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)])

    # \chi^3\left(\vec{r}\right)
    omega = chi3

    # \frac{\nabla \varepsilon_r\left(\vec{r}\right) / \varepsilon_r\left(\vec{r}\right)}
    etax = Data.Array([ 0*exp(-(y_g(iy, dy, ux[:, :, :]) - ly / 2)^2 - (z_g(iz, dz, ux[:, :, :]) - lz / 2)^2) / epsilon[ix,iy,iz] for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)])
    etay = Data.Array([ -2*(y_g(iy, dy, ux[:, :, :]) - ly / 2)*exp(-(y_g(iy, dy, ux[:, :, :]) - ly / 2)^2 - (z_g(iz, dz, ux[:, :, :]) - lz / 2)^2) / epsilon[ix,iy,iz] for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)])
    etaz = Data.Array([ -2*(z_g(iz, dz, ux[:, :, :]) - lz / 2)*exp(-(y_g(iy, dy, ux[:, :, :]) - ly / 2)^2 - (z_g(iz, dz, ux[:, :, :]) - lz / 2)^2) / epsilon[ix,iy,iz] for ix = 1:size(ux, 1), iy = 1:size(ux, 2), iz = 1:size(ux, 3)])
    
    # laod epsilon to gpu for alpha update
    epsilon = Data.Array(epsilon)
    
    if do_visu
        ENV["GKSwstype"] = "nul"
        if (me == 0)
            if isdir(example_name) == false
                mkdir(example_name)
            end
            loadpath = example_name * "/"
            anim = Animation(loadpath, String[])
            println("Animation directory: $(anim.dir)")
        end
        nx_v, ny_v, nz_v = (nx - 2) * dims[1], (ny - 2) * dims[2], (nz - 2) * dims[3]
        if (4 * nx_v * ny_v * nz_v * sizeof(Data.Number) > 0.8 * Sys.free_memory())
            error("Not enough memory for visualization.")
        end
        u_v = zeros(3, nx_v, ny_v, nz_v)
        u_tmp = zeros(nx_v, ny_v, nz_v)
        ux_inn = zeros(nx - 2, ny - 2, nz - 2)
        uy_inn = zeros(nx - 2, ny - 2, nz - 2)
        uz_inn = zeros(nx - 2, ny - 2, nz - 2)
        xi_g, yi_g, zi_g = LinRange(0, lx, nx_v), LinRange(0, ly, ny_v), LinRange(0, lz, nz_v) # inner points only
        iframe = 0
    end

    for it = 1:nt
        # non static coefficient function  
        # call compute alpha kernel
        @parallel (1:size(ux, 1), 1:size(ux, 2), 1:size(ux, 3)) compute_alpha!(alpha, ux, uy, uz, _mu0, _epsilon0, chi3, epsilon)

        # update other coefficients with array programming
        gamma   .= gamma .* alpha
        beta    .= beta .* alpha 

        # update v with the right kernels
        @parallel (1:size(ux, 1), 1:size(ux, 2), 1:size(ux, 3)) update_vecv_nabla2!(ux, uy, uz, vx, vy, vz, dt, _dx2, _dy2, _dz2, alpha)
        @parallel (1:size(ux, 1)-2, 1:size(ux, 2)-2, 1:size(ux, 3)-2) update_vecv_sigma!(ux, uy, uz, vx, vy, vz, dt, beta, gamma)
        @parallel (1:size(ux, 1), 1:size(ux, 2), 1:size(ux, 3)) update_vecv_varepsilon!(ux, uy, uz, vx, vy, vz, dt, _dx_2, _dy_2, _dz_2, alpha, etax, etay, etaz)
        @parallel (1:size(ux, 1), 1:size(ux, 2), 1:size(ux, 3)) update_vecv_chi3!(ux, uy, uz, vx, vy, vz, dt, _dx, _dy, _dz, _dx_2, _dy_2, _dz_2, alpha, delta, omega)
        
        # update u
        @hide_communication b_width begin
            @parallel update_vecu!(ux, uy, uz, vx, vy, vz, dt)
            update_halo!(ux, uy, uz)
        end

        if do_visu && (it % nvis == 0)
            ux_inn .= Array(ux)[2:end-1, 2:end-1, 2:end-1]
            uy_inn .= Array(uy)[2:end-1, 2:end-1, 2:end-1]
            uz_inn .= Array(uz)[2:end-1, 2:end-1, 2:end-1]
            gather!(ux_inn, u_tmp)
            u_v[1, :, :, :] .= u_tmp
            gather!(uy_inn, u_tmp)
            u_v[2, :, :, :] .= u_tmp
            gather!(uz_inn, u_tmp)
            u_v[3, :, :, :] .= u_tmp

            if me == 0
                plt = heatmap(xi_g, yi_g, u_v[2, :, :, ceil(Int, nz_g() / 2)]'; xlims=(xi_g[1], xi_g[end]), ylims=(yi_g[1], yi_g[end]), clims=(-1, 1), aspect_ratio=1, c=:turbo)
                png(plt, example_name * @sprintf("/%05d.png", iframe += 1))
                save_array(example_name * @sprintf("/out_T_%05d", iframe), convert.(Float32, u_v))
            end
        end
    end
    finalize_global_grid()
    return nothing
end


nonlinear_vecE_3D(do_visu=true)