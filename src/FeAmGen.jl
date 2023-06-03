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
using OrderedCollections
using Pipe
using Pkg
using PyCall
using SHA
using SymEngine
using YAML

export digest_seed_proc
export generate_amp 
export canonicalize_amp
export construct_den_topology

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
include("Topology.jl")
include("Visual.jl")


###################
function __init__()
###################
  return nothing
end # function __init__


end # module FeAmGen
