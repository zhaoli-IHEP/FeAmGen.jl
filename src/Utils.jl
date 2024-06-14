# Copyright (c) 2023 Fenyutanchan <fenyutanchan@vip.qq.com>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

#################################################
"""
    find_fermion_loops( tex_file::String )::Vector{Vector{String}}

Give the particle names of the fermion loop in the Feynman diagram.
"""
function find_fermion_loops( tex_file::String )::Vector{Vector{String}}
#################################################
  @assert isfile(tex_file) "File not found: $tex_file."
  v_line_list = filter!( !isempty, readlines( tex_file ) )
  for line_index ∈ eachindex(v_line_list)
    while startswith( v_line_list[line_index], r"\s")
      v_line_list[line_index] = v_line_list[line_index][begin+1:end]
    end # while
    while endswith( v_line_list[line_index], r"\s")
      v_line_list[line_index] = v_line_list[line_index][begin:end-1]
    end # while
  end # for line_index
  v_line_list = filter!( startswith(r"v[1-9]\d*"), v_line_list )
  v_line_list = filter!( endswith(r"v[1-9]\d*,"), v_line_list )
  v_line_list = filter!( contains("fermion"), v_line_list )
  isempty(v_line_list) && return [String[]]

  v_pair_list = Vector{Int}[]
  particle_name_list = String[]
  for v_line ∈ v_line_list
    v_pair = map( v_str_range->parse(Int,match(r"\d+",v_line[v_str_range]).match), findall(r"v[1-9]\d*",v_line) )
    particle_name = match( r"edge label' = \\\(\w+?\\\),", v_line ).match
    particle_name = replace( particle_name, "edge label' = \\(" => "", "\\)," => "" )
    push!( v_pair_list, v_pair )
    push!( particle_name_list, particle_name )
  end # for v_line

  while true
    all_v_id_list = vcat( v_pair_list... )
    v_id_list = (sort!∘unique)( all_v_id_list )
    v_count_list = map( v_id->count(==(v_id),all_v_id_list), v_id_list )
    @assert (iszero∘count)(≥(3), v_count_list) "More than two fermion lines adjacent to one vertex: $tex_file."
    end_v_id_indices = findall( ≠(2), v_count_list )
    isempty(end_v_id_indices) && break

    end_v_id = v_id_list[end_v_id_indices]
    to_be_deleted_indices = findall( v_pair->!isempty(end_v_id∩v_pair), v_pair_list )
    deleteat!( v_pair_list, to_be_deleted_indices )
    deleteat!( particle_name_list, to_be_deleted_indices )
  end # while

  isempty(v_pair_list) && return [String[]]

  loop_list = Vector{String}[]
  while !isempty(v_pair_list)
    v_pair = popfirst!(v_pair_list)
    loop = [popfirst!(particle_name_list)]

    while true
      selected_index = findfirst( ==(last(v_pair)) ∘ first, v_pair_list )
      isnothing(selected_index) && break
      push!( v_pair, (last ∘ popat!)( v_pair_list, selected_index ) )
      push!( loop, popat!( particle_name_list, selected_index ) )
    end # while

    push!( loop_list, loop )
  end # while

  return loop_list
end # function find_fermion_loops



#################################################
"""
    gen_shifted_amp( top_file::String; path_type::Symbol=:dir )::Vector{String}

Generate the shifted amplitudes according to the topology file.
"""
function gen_shifted_amp( top_file::String; path_type::Symbol=:dir )::Vector{String}
#################################################
  @assert ispath( top_file ) "Path not found: $top_file."

  shifted_amp_dir = "_shifted_amps"
  if path_type == :dir
    @assert isdir( top_file ) "Directory not found: $top_file."
    shifted_amp_dir = top_file * shifted_amp_dir
    bk_mkdir( shifted_amp_dir )
    top_file_list = filter!( endswith(".jld2"), readdir( top_file; join=true) )
    file_list = map( top_file -> gen_shifted_amp(top_file; path_type=:file), top_file_list )
    return union( file_list... )
  end # if
  @assert isfile( top_file ) "File not found: $top_file."
  @assert endswith( top_file, ".jld2" ) "Invalid file extension: $top_file."
  top_info = load( top_file )
  @assert (isempty ∘ symdiff)( keys(top_info), ["kinematic_relation", "covering_amplitudes", "external_momenta", "loop_momenta", "denominators"] ) "Invalid topology file: $top_file."

  top_dir = (dirname ∘ abspath)( top_file )
  shifted_amp_dir = top_dir * shifted_amp_dir
  if !isdir( shifted_amp_dir )
    rm( shifted_amp_dir; recursive=true, force=true )
    mkdir( shifted_amp_dir )
  end # if

  shifted_amp_list = String[]
  for (amp_file, mom_shift_dict) ∈ top_info["covering_amplitudes"]
    original_amp_file = (abspath ∘ joinpath)( top_dir, amp_file )
    @info "Input amplitude: $original_amp_file"
    @assert isfile( original_amp_file ) "File not found: $original_amp_file."
    amp_info = load( original_amp_file )
    target_amp_file = joinpath( shifted_amp_dir, basename( original_amp_file )[begin:end-5] * "_" * basename( top_file ) )

    mom_shift = Dict{Basic, Basic}()
    for (k, v) ∈ mom_shift_dict
      mom_shift[Basic(k)] = Basic(v)
    end # for (k, v)

    loop_den_list = map(Basic, amp_info["loop_den_list"])
    amp_lorentz_list = map(Basic, amp_info["amp_lorentz_list"])

    for (ii, loop_den) ∈ enumerate(loop_den_list)
      loop_den_list[ii] = subs( loop_den, mom_shift )
    end # for (ii, loop_den)
    loop_den_list = normalize_loop_mom( loop_den_list )
    @assert (isempty∘setdiff)(loop_den_list, to_Basic(top_info["denominators"])) begin
      """
      The loop denominators in the amplitude file are not consistent with the topology file.
        shifted_loop_dens: $loop_den_list
        topology_dens: $(top_info["denominators"])
      """
    end # @assert
    for (ii, amp_lorentz) ∈ enumerate(amp_lorentz_list)
      amp_lorentz_list[ii] = subs( amp_lorentz, mom_shift )
    end # for (ii, amp_lorent)

    amp_info["loop_den_list"] = map(string, loop_den_list)
    amp_info["amp_lorentz_list"] = map(string, amp_lorentz_list)
    # amp_info["corresponding_topology"] = relative_path( target_amp_file, top_file )

    jldopen( target_amp_file, "w" ) do io
      for (k, v) ∈ amp_info
        io[k] = v
      end # for (k, v)
      io["topology"] = top_info
    end # jldopen
    
    push!( shifted_amp_list, target_amp_file )
    @info "Shifted amplitudes generated: $target_amp_file"
  end # for (amp_file, mom_shift)

  return shifted_amp_list
end # function gen_shifted_amp



remove_recover_underscore( expr::Basic ) = remove_recover_underscore( [expr] )
function remove_recover_underscore(
  exprs::Vector{Basic}
)::Tuple{
  Dict{String,String}, # removed_underscore_dict
  Dict{String,String} # recover_underscore_dict
}
  with_underscore_symbols = map( string, free_symbols( exprs ) )
  filter!( contains("_"), with_underscore_symbols )
  removed_underscore_symbols = map( x -> replace(x,"_"=>"") , with_underscore_symbols )
  unique_removed_underscore_symbols = unique( removed_underscore_symbols )
  for unique_removed_underscore_symbol ∈ unique_removed_underscore_symbols
    positions = findall( ==(unique_removed_underscore_symbol), removed_underscore_symbols )
    length(positions) == 1 && continue
    for (ii, position) ∈ enumerate(positions)
      removed_underscore_symbols[position] *= "repeated$ii"
    end # for
  end # for unique_removed_underscore_symbol

  removed_underscore_dict = Dict( with_underscore_symbols .=> removed_underscore_symbols )
  recovered_underscore_dict = Dict( removed_underscore_symbols .=> with_underscore_symbols )

  return removed_underscore_dict, recovered_underscore_dict
end # function remove_recover_underscore
