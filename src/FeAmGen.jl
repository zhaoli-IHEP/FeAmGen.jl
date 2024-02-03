__precompile__()

"""
    module FeAmGen

FeAmGen represents `Fe`ynman diagram and `Am`plitude `Gen`erator.
"""
module FeAmGen

using AbstractAlgebra
using FeynUtils
using FeynmanDenominators
using Combinatorics
using Dates
using JLD2
using nauty_jll
using OrderedCollections
using PakAlgorithm
using Pipe
using Pkg
using PyCall
using SymEngine
using YAML

import FORM_jll

# Main APIs
export digest_seed_proc
export generate_amp 
export construct_den_topology, construct_topology

# Extra important functionals
export canonicalize_amp
export find_fermion_loops
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
include("External.jl")
include("Seed.jl")
include("SimpleDigest.jl")
include("ToMathematicaForm.jl")
include("Topology.jl")
include("Visual.jl")
include("Utils.jl")


###################
function __init__()
###################
  return nothing
end # function __init__


end # module FeAmGen
