# push!(LOAD_PATH,"../src/")
include("make_general_1.jl")
include("make_external_1.jl")
include("make_general_2.jl")

# deploydocs(
#     deps = Deps.pip("pygments", "mkdocs", "python-markdown-math")
# )

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
