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
  isempty(v_line_list) && return Int[]

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

    end_v_id = v_id_list[v_id_list]
    to_be_deleted_indices = findall( v_pair->end_v_id∈v_pair, v_pair_list )
    deleteat!( v_pair_list, to_be_deleted_indices )
    deleteat!( particle_name_list, to_be_deleted_indices )
  end # while

  isempty(v_pair_list) && return [String[]]
  sort!(v_pair_list; by=first)

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
