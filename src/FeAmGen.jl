__precompile__()

"""
    module FeAmGen

FeAmGen represents `Fe`ynman diagram and `Am`plitude `Gen`erator.
"""
module FeAmGen

using AbstractAlgebra
using FeynUtils
using Combinatorics
using Dates
using FORM_jll
using JLD2
using nauty_jll
using OrderedCollections
using Pipe
using Pkg
using PyCall
using SHA
using SymEngine
using YAML

# Main APIs
export construct_den_topology
export digest_seed_proc
export generate_amp 

# Extra important functionals
export canonicalize_amp
export is_planar
export to_m_file

include("Universe.jl")
include("Graph.jl")
include("Canon.jl")
include("Digest.jl")
include("FeynmanDiagram.jl")
include("FORMS.jl")
include("GenAmp.jl")
include("Kin.jl")
include("QGRAF.jl")
include("Seed.jl")
include("SimpleDigest.jl")
include("ToMathematicaForm.jl")
include("Topology.jl")
include("Visual.jl")


###################
function __init__()
###################
  return nothing
end # function __init__


end # module FeAmGen
