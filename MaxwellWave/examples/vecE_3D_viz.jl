# instantiate for running on daint
using Documenter
using Printf
using GLMakie
using LinearAlgebra


include("../src/auxiliary.jl")

@doc raw"""
    vizme()

Generates a 3D gifs and plots from the data generated by any vec 3D exampel files.
The settings in the function have to by changed to match the choosen script for a functioning script.
"""
@views function vizme()

    save_folder = "docs/" 
    # this had to be changed to chose the script
    example_name = "homogeneous_dirichlet_vecE_3D"

    # include the physics and numerics settings for plotting
    include( example_name * "_settings.jl")

    # number of ranks from mpi run, this had to be changed to match the mpi config
    dims = [2,2,2]

    # sparsitiy for vector plot
    sparse = [6,4,4]

    # numerics, this had to be changed to match the script 
    nx_v,ny_v,nz_v = (nx-2)*dims[1],(ny-2)*dims[2],(nz-2)*dims[3]
    u       = zeros(Float32,3,nx_v,ny_v,nz_v)
    xc      = LinRange(0, lx, nx_v)
    yc      = LinRange(0, ly, ny_v)
    zc      = LinRange(0, lz, nz_v)
    dx      = lx/(nx_v-1)
    dy      = ly/(ny_v-1)
    dz      = lz/(nz_v-1)
    dt = min(dx, dy, dy) / sqrt(c2) / 8
    
    # physical constants
    # permitivity, this had to be changed to match the script
    epsilon = [(1 + exp(-(yc[iy] - ly / 2)^2 - (zc[iz] - lz / 2)^2)) for ix = 1:size(u, 2), iy = 1:size(u, 3), iz = 1:size(u, 4)]
    # conductivity, this had to be changed to match the script
    sigma = [(1 - exp(-0.5*(yc[iy] - ly / 2)^2 - 0.5*(zc[iz] - lz / 2)^2)) for ix = 1:size(u, 2), iy = 1:size(u, 3), iz = 1:size(u, 4)]
    
    # chi3, this had to be changed to match the script 
    chi3 = [0.01*exp(-(xc[ix] - lx / 2)^2 -(yc[iy] - ly / 2)^2 - (zc[iz] - lz / 2)^2) for ix = 1:size(u, 2), iy = 1:size(u, 3), iz = 1:size(u, 4)]

    # initial field to plot, this had to be changed to match the script
    u_pulse_shape = [2 * sech(0.5 * (k0 * (xc[ix] - lx / 2))) * exp(-(yc[iy] - ly / 2)^2 / (2 * sigma2[2]) - (zc[iz] - lz / 2)^2 / (2 * sigma2[3])) for ix = 1:size(u, 2), iy = 1:size(u, 3), iz = 1:size(u, 4)]
    ux = p0[1] .* u_pulse_shape
    uy = p0[2] .* u_pulse_shape
    uz = p0[3] .* u_pulse_shape   

    v_pulse_shape = [(2 * sech(0.5 * (k0 * (xc[ix] - lx / 2) - w0 * dt / 2)) - 2 * sech(0.5 * (k0 * (xc[ix] - lx / 2) + w0 * dt / 2))) / dt * exp(-(yc[iy] - ly / 2)^2 / (2 * sigma2[2]) - (zc[iz] - lz / 2)^2 / (2 * sigma2[3])) for ix = 1:size(u, 2), iy = 1:size(u, 3), iz = 1:size(u, 4)]
    vx = p0[1] .* v_pulse_shape
    vy = p0[2] .* v_pulse_shape
    vz = p0[3] .* v_pulse_shape




    
    # create figure 
    fig     = Figure(resolution=(1600,1000),fontsize=24,backgroundcolor = "gray")


    # create vector quiver plot gif
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Electric Field",xlabel="lx",ylabel="ly",zlabel="lz")
    load_array(example_name * "/out_T_00001",u)

    
    ps = [Point3f(x, y, z) for x in xc[1:sparse[1]:end] for y in yc[1:sparse[2]:end] for z in zc[1:sparse[3]:end]]
    ns = map(p -> Vec3f(u[1,Int(round(p[1]/dx+1)),Int(round(p[2]/dy+1)),Int(round(p[3]/dz+1))], u[2,Int(round(p[1]/dx+1)),Int(round(p[2]/dy+1)),Int(round(p[3]/dz+1))], u[3,Int(round(p[1]/dx+1)),Int(round(p[2]/dy+1)),Int(round(p[3]/dz+1))]), ps)
    lengths = norm.(ns)
    plt = arrows!(
        ps, ns, fxaa=true, # turn on anti-aliasing
        color=lengths,
        normalize = false,
        colormap = :turbo,  
        arrowsize = Vec3f(0.1, 0.1, 0.02),
        linewidth = 0.08
    )

    record(fig, save_folder * example_name * "_vecplt.gif", 2:nt/nvis; framerate = 4) do iframe
        empty!(fig)

        load_array(example_name * @sprintf("/out_T_%05d",iframe),u)
        ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Electric Field",xlabel="lx",ylabel="ly",zlabel="lz")
        ps = [Point3f(x, y, z) for x in xc[1:sparse[1]:end] for y in yc[1:sparse[2]:end] for z in zc[1:sparse[3]:end]]
        ns = map(p -> Vec3f(u[1,Int(round(p[1]/dx+1)),Int(round(p[2]/dy+1)),Int(round(p[3]/dz+1))], u[2,Int(round(p[1]/dx+1)),Int(round(p[2]/dy+1)),Int(round(p[3]/dz+1))], u[3,Int(round(p[1]/dx+1)),Int(round(p[2]/dy+1)),Int(round(p[3]/dz+1))]), ps)
        lengths = norm.(ns)
        plt = arrows!(
            ps, ns, fxaa=true, # turn on anti-aliasing
            color=lengths,
            normalize = false,
            colormap = :turbo,  
            arrowsize = Vec3f(0.1, 0.1, 0.02),
            linewidth = 0.08
        )
    end



    # create Ex slice plot gif
    empty!(fig)
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Electric Field x-Component",xlabel="lx",ylabel="ly",zlabel="lz")
    load_array(example_name * "/out_T_00001",u)
    plt = volumeslices!(ax, xc, yc, zc, u[1,:,:,:],
                        colormap = :turbo,  
                        colorrange=(-1, 1),
                        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))
    record(fig, save_folder * example_name * "_slicepltx.gif", 2:nt/nvis; framerate = 4) do iframe
        empty!(fig)

        load_array(example_name * @sprintf("/out_T_%05d",iframe),u)
        ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Electric Field x-Component",xlabel="lx",ylabel="ly",zlabel="lz")
        plt = volumeslices!(ax, xc, yc, zc, u[1,:,:,:],
                            colormap = :turbo,  
                            colorrange=(-1, 1),
                            transparency = true
                            )
        plt[:update_yz][](ceil(Int64,nx_v/2))
        plt[:update_xz][](ceil(Int64,ny_v/2))
        plt[:update_xy][](ceil(Int64,nz_v/2))
    end



    # create Ey slice plot gif
    empty!(fig)
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Electric Field y-Component",xlabel="lx",ylabel="ly",zlabel="lz")
    load_array(example_name * "/out_T_00001",u)
    plt = volumeslices!(ax, xc, yc, zc, u[2,:,:,:],
                        colormap = :turbo,  
                        colorrange=(-1, 1),
                        transparency = true
    )
    
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))
    record(fig, save_folder * example_name * "_sliceplty.gif", 2:nt/nvis; framerate = 4) do iframe
        empty!(fig)

        load_array(example_name * @sprintf("/out_T_%05d",iframe),u)
        ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Electric Field y-Component",xlabel="lx",ylabel="ly",zlabel="lz")
        plt = volumeslices!(ax, xc, yc, zc, u[2,:,:,:],
                            colormap = :turbo,  
                            colorrange=(-1, 1),
                            transparency = true
                            )
        plt[:update_yz][](ceil(Int64,nx_v/2))
        plt[:update_xz][](ceil(Int64,ny_v/2))
        plt[:update_xy][](ceil(Int64,nz_v/2))
    end



    # create Ez slice plot gif
    empty!(fig)
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Electric Field z-Component",xlabel="lx",ylabel="ly",zlabel="lz")
    load_array(example_name * "/out_T_00001",u)
    plt = volumeslices!(ax, xc, yc, zc, u[3,:,:,:],
                        colormap = :turbo,  
                        colorrange=(-1, 1),
                        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))
    record(fig, save_folder * example_name * "_slicepltz.gif", 2:nt/nvis; framerate = 4) do iframe
        empty!(fig)

        load_array(example_name * @sprintf("/out_T_%05d",iframe),u)
        ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Electric Field z-Component",xlabel="lx",ylabel="ly",zlabel="lz")
        plt = volumeslices!(ax, xc, yc, zc, u[3,:,:,:],
                            colormap = :turbo,  
                            colorrange=(-1, 1),
                            transparency = true
                            )
        plt[:update_yz][](ceil(Int64,nx_v/2))
        plt[:update_xz][](ceil(Int64,ny_v/2))
        plt[:update_xy][](ceil(Int64,nz_v/2))
    end



    empty!(fig)
    # material permitivity plot
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Permeability",xlabel="lx",ylabel="ly",zlabel="lz")

    plt = volumeslices!(ax, xc, yc, zc, epsilon,
        colormap = :turbo,  
        colorrange=(0, 2*epsilon0),
        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))

    save(save_folder * example_name * "_epsilon.png", fig)



    empty!(fig)
    # material conductivity plot
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Conductivity",xlabel="lx",ylabel="ly",zlabel="lz")

    plt = volumeslices!(ax, xc, yc, zc, sigma,
        colormap = :turbo,  
        colorrange=(0, 2*epsilon0),
        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))

    save(save_folder * example_name * "_sigma.png", fig)


    empty!(fig)
    # material conductivity plot
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Chi3",xlabel="lx",ylabel="ly",zlabel="lz")

    plt = volumeslices!(ax, xc, yc, zc, chi3,
        colormap = :turbo,  
        colorrange=(0, 0.01),
        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))

    save(save_folder * example_name * "_chi3.png", fig)


    empty!(fig)
    # material conductivity plot
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Initial Field x",xlabel="lx",ylabel="ly",zlabel="lz")

    plt = volumeslices!(ax, xc, yc, zc, ux,
        colormap = :turbo,  
        colorrange=(-1, 1),
        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))

    save(save_folder * example_name * "_initux.png", fig)



    empty!(fig)
    # material conductivity plot
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Initial Field y",xlabel="lx",ylabel="ly",zlabel="lz")

    plt = volumeslices!(ax, xc, yc, zc, uy,
        colormap = :turbo,  
        colorrange=(-1, 1),
        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))

    save(save_folder * example_name * "_inituy.png", fig)



    empty!(fig)
    # material conductivity plot
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Initial Field z",xlabel="lx",ylabel="ly",zlabel="lz")

    plt = volumeslices!(ax, xc, yc, zc, uz,
        colormap = :turbo,  
        colorrange=(-1, 1),
        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))

    save(save_folder * example_name * "_inituz.png", fig)

    

    empty!(fig)
    # material conductivity plot
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Initial Velocity x",xlabel="lx",ylabel="ly",zlabel="lz")

    plt = volumeslices!(ax, xc, yc, zc, vx,
        colormap = :turbo,  
        colorrange=(-1, 1),
        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))

    save(save_folder * example_name * "_initvx.png", fig)

    
    
    empty!(fig)
    # material conductivity plot
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Initial Velocity y",xlabel="lx",ylabel="ly",zlabel="lz")

    plt = volumeslices!(ax, xc, yc, zc, vy,
        colormap = :turbo,  
        colorrange=(-1, 1),
        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))

    save(save_folder * example_name * "_initvy.png", fig)

    
    empty!(fig)
    # material conductivity plot
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Initial Velocity z",xlabel="lx",ylabel="ly",zlabel="lz")

    plt = volumeslices!(ax, xc, yc, zc, vz,
        colormap = :turbo,  
        colorrange=(-1, 1),
        transparency = true
    )
    plt[:update_yz][](ceil(Int64,nx_v/2))
    plt[:update_xz][](ceil(Int64,ny_v/2))
    plt[:update_xy][](ceil(Int64,nz_v/2))

    save(save_folder * example_name * "_initvz.png", fig)

    
    return nothing
end


vizme()