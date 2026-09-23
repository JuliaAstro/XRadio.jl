using Documenter
using XRadio

DocMeta.setdocmeta!(XRadio, :DocTestSetup, :(using XRadio); recursive = true)

makedocs(;
    sitename = "XRadio.jl",
    modules = [XRadio],
    checkdocs = :public,
    pages = [
        "Home" => "index.md",
        "Selecting data" => "selection.md",
        "Modifying data" => "modifying.md",
        "Writing stores" => "writing.md",
        "Building data" => "building.md",
        "Coming from xradio" => "xradio.md",
        "Reference" => "reference.md",
    ],
)
