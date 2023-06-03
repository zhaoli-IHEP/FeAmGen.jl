push!(LOAD_PATH,"../src/")

using Documenter, FeAmGen

makedocs( 
  modules = [FeAmGen],
  format = Documenter.HTML(
        # Disable pretty URLs during manual testing
        prettyurls = get(ENV, "CI", nothing) == "true",
  ),
  sitename = "FeAmGen.jl Documentation", 
  pages = Any[
    "Home" => "index.md",
    "Types" => "types.md",
    "Front-end functions" => "front_end.md",
    "Functions" => [ 
       "Functions in Message.jl" => "FunctionsMessage.md",
       "Functions in Extra.jl" => "FunctionsExtra.md",
       "Functions in Universe.jl" => "FunctionsUniverse.md",
       "Functions in SimpleDigest.jl" => "FunctionsSimpleDigest.md",
       "Functions in Digest.jl" => "FunctionsDigest.md",
       "Functions in Kin.jl" => "FunctionsKin.md",
       "Functions in FORMS.jl" => "FunctionsFORMS.md",
       "Functions in Visual.jl" => "FunctionsVisual.md",
       "Functions in Converter.jl" => "FunctionsConverter.md",
       "Functions in FeynmanDiagram.jl" => "FunctionsFeynmanDiagram.md"
    ]
  ]
)

