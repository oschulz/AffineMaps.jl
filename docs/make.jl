# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using AffineTransformations

# Doctest setup
DocMeta.setdocmeta!(
    AffineTransformations,
    :DocTestSetup,
    :(using AffineTransformations);
    recursive=true,
)

makedocs(
    sitename = "AffineTransformations",
    modules = [AffineTransformations],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://oschulz.github.io/AffineTransformations.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    strict = !("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/oschulz/AffineTransformations.jl.git",
    forcepush = true,
    push_preview = true,
)
