using Documenter, MaxwellWave

makedocs(
    modules = [MaxwellWave],
    format = Documenter.HTML(
        prettyurls = false,
        canonical = "http://MaxwellWave.almaeder.com/stable/",
    ),
    sitename = "MaxwellWave.jl",
    pages = [
        "index.md",
        "wave2D.md",
        "wave3D.md",
        "wave3DnonlinearChi3.md",
        "auxiliary.md",
    ],
)
deploydocs(
    repo = "github.com/lheimabch/MaederHeimbach/MaxwellWave.jl.git",devbranch = "heimbach"
)