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
    ],
)