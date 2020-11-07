push!(LOAD_PATH,"../src/")
using Documenter
using OceanAcoustics

makedocs(
    sitename = "Ocean Acoustics",
    format = Documenter.HTML(),
    modules = [OceanAcoustics]
)

# deploydocs(
#     deps = Deps.pip("pygments", "mkdocs", "python-markdown-math")
# )

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
