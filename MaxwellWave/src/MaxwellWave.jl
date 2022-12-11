module MaxwellWave

module Auxiliary
include("../src/auxiliary.jl")
end

module Wave2D
using ParallelStencil, ParallelStencil.FiniteDifferences3D
@init_parallel_stencil(Threads, Float64, 3)
include("../src/wave2D.jl")
end

module Wave3D
using ParallelStencil, ParallelStencil.FiniteDifferences3D
@init_parallel_stencil(Threads, Float64, 3)
include("../src/wave3D.jl")
end

end # module MaxwellWave
