using Documenter

using Plots, ParallelStencil, ParallelStencil.FiniteDifferences2D, Printf 
using ImplicitGlobalGrid
# import MPI

# decide if one uses the gpu or cpu
const USE_GPU = false
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end


include("../src/wave2D.jl")
include("../src/auxiliary.jl")

@doc raw"""
    homogeneous_dirichlet_Ez_2D(; do_visu=false)

Solves the two dimensional wave equation problem:

```math
\nabla^2 E_z - \frac{1}{c^2} \frac{\partial^2 }{\partial t^2} E_z = 0, \; E_z \in \Omega
```
```math
E_z = 0, \; E_z \in \partial \Omega
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
function homogeneous_dirichlet_Ez_2D(; do_visu=false)
    # physical parameters 
    lx = 10
    ly = lx
    c2 = 1
    lambda = 1
    k = 2 * pi / lambda
    w = k * sqrt(c2)

    # numerical parameters
    nx = 100
    ny = nx
    nt = 2000
    nvis = 50

    # initailize glo bal grid and such MPI
    me, dims = init_global_grid(nx, ny, 1)
    # optimization for communication
    b_width = (8, 8, 4)

    # derived numerics
    dx = lx / nx_g()
    dy = ly / ny_g()
    _dx = 1 / dx
    _dy = 1 / dy
    dt = min(dx, dy) / sqrt(c2) / 2

    # array initialization 
    u = @zeros(nx, ny)
    u .= Data.Array([exp(-(x_g(ix, dx, u) - lx / 2)^2 - (y_g(iy, dy, u) - ly / 2)^2) * cos(k * x_g(ix, dx, u)) for ix = 1:size(u, 1), iy = 1:size(u, 2)])

    update_halo!(u)
    v = @zeros(nx, ny)
    v .= Data.Array([-exp(-(x_g(ix, dx, u) - lx / 2)^2 - (y_g(iy, dy, u) - ly / 2)^2) * w * sin(k * x_g(ix, dx, u)) for ix = 1:size(u, 1), iy = 1:size(u, 2)])


    if do_visu
        ENV["GKSwstype"]="nul"
        if (me==0) if isdir("homogeneous_dirichlet_Ez_out")==false mkdir("homogeneous_dirichlet_Ez_out") end; loadpath="homogeneous_dirichlet_Ez_out/"; anim=Animation(loadpath,String[]); println("Animation directory: $(anim.dir)") end
        nx_v, ny_v = (nx - 2) * dims[1], (ny - 2) * dims[2]
        if (nx_v * ny_v * sizeof(Data.Number) > 0.8 * Sys.free_memory())
            error("Not enough memory for visualization.")
        end
        u_v = zeros(nx_v, ny_v)
        u_inn = zeros(nx - 2, ny - 2)
        xi_g, yi_g = LinRange(0, lx, nx_v), LinRange(0, ly, ny_v) # inner points only
        iframe = 0
    end

    for it = 1:nt
        @parallel (1:size(u, 1), 1:size(u, 2)) update_v!(u, v, dt, _dx, _dy, c2)
        @hide_communication b_width begin
            @parallel (1:size(u, 1) - 2, 1:size(u, 2) - 2) update_u!(u, v, dt)
            update_halo!(u)
        end
        if do_visu && (it % nvis == 0)
            u_inn .= Array(u)[2:end-1, 2:end-1]
            gather!(u_inn, u_v)
            if me == 0
                plt = heatmap(xi_g, yi_g, u_v'; xlims=(xi_g[1], xi_g[end]), ylims=(yi_g[1], yi_g[end]), clims=(-1, 1), aspect_ratio=1, c=:turbo)
                png(plt,@sprintf("homogeneous_dirichlet_Ez_out/%04d.png",iframe+=1))
                save_array(@sprintf("homogeneous_dirichlet_Ez_out/out_T_%04d",iframe),convert.(Float32,u_v))
            end
        end
    end


    finalize_global_grid()
    return nothing
end

homogeneous_dirichlet_Ez_2D(do_visu=true)