using ParallelStencil, ParallelStencil.FiniteDifferences3D
using Test
@testset "wave2d" begin
    include("../src/wave2D.jl")
    @init_parallel_stencil(Threads, Float64, 3)

    u_init = rand(5, 5)
    u_gold = copy(u_init)
    u_init = Data.Array(u_init)
    v_init = rand(5, 5)
    v_gold = copy(v_init)
    v_init = Data.Array(v_init)
    dt = rand(1)[1]
    _dx = rand(1)[1]
    _dy = rand(1)[1]
    c2 = rand(1)[1] 
    @parallel (1:size(u_init, 1)-2, 1:size(u_init, 2)-2) update_u!(u_init, v_init, dt)
    @parallel (1:size(u_init, 1), 1:size(u_init, 2)) update_v!(u_init, v_init, dt, _dx, _dy, c2)
    u_final = Array(u_init)
    v_final = Array(v_init)
    @test all(u_gold[1:end,1] .== u_final[1:end,1])
    @test all(u_gold[1:end,end] .== u_final[1:end,end])
    @test all(u_gold[1,1:end] .== u_final[1,1:end])
    @test all(u_gold[end,1:end] .== u_final[end,1:end])
    @test all(v_gold[1:end,1] .== v_final[1:end,1])
    @test all(v_gold[1:end,end] .== v_final[1:end,end])
    @test all(v_gold[1,1:end] .== v_final[1,1:end])
    @test all(v_gold[end,1:end] .== v_final[end,1:end])

    u_init = LinRange(0,rand(1)[1]+1,5)*LinRange(0,rand(1)[1]+1,5)'
    u_init = Data.Array(u_init)
    v_init = rand(5, 5)
    v_gold = copy(v_init)
    v_init = Data.Array(v_init)
    dt = rand(1)[1]
    _dx = rand(1)[1]
    _dy = rand(1)[1]
    c2 = rand(1)[1] 
    @parallel (1:size(u_init, 1), 1:size(u_init, 2)) update_v!(u_init, v_init, dt, _dx, _dy, c2)
    v_final = Array(v_init)
    @test all(v_final .== v_gold)

    u_init = rand(5, 5)
    u_gold = copy(u_init)
    u_init = Data.Array(u_init)
    v_init = zeros(5,5)
    v_init = Data.Array(v_init)
    dt = rand(1)[1]
    @parallel (1:size(u_init, 1)-2, 1:size(u_init, 2)-2) update_u!(u_init, v_init, dt)
    u_final = Array(u_init)
    @test all(u_final .== u_gold)
end
@testset "wave3D" begin
    @test all(1 .== 1)
end
@testset "auxiliary" begin
    include("../src/auxiliary.jl")
    A = rand(3,3)
    B = similar(A)
    save_array("bin_io",A)
    load_array("bin_io",B)
    @test all(A .== B)
end