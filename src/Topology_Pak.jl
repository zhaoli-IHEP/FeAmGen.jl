function construct_den_topology(
  ::Val{:PakAlgorithm},
  den_collection_list::Vector{FeynmanDenominatorCollection},
  kinematic_relations::Dict{Basic, Basic},
  find_external_momentum_shifts::Bool
)::Tuple{
  Vector{FeynmanDenominatorCollection},
  Vector{Dict{Int, Dict{Basic, Basic}}}
}
  result_den_collection_list = FeynmanDenominatorCollection[]
  result_repl_rules_list = Dict{Int, Dict{Basic, Basic}}[]
  remaining_indices = (collect ∘ eachindex)(den_collection_list)

  @info "Minimalizing the topologies with Pak algorithm."
  while !isempty(remaining_indices)
    @info "Number of remaining topologies: $(length(remaining_indices))."

    _, index_index = findmax(length, den_collection_list[remaining_indices])
    largest_index = remaining_indices[index_index]
    # index_index, largest_index = 1, first(remaining_indices)
    while true
      tmp_den_collection = den_collection_list[largest_index]
      tmp_index_index = findnext(
        index -> tmp_den_collection ≠ den_collection_list[index] &&
          !is_Pak_equivalent(
            tmp_den_collection, den_collection_list[index];
            SP_replacements=kinematic_relations
          ) && den_collection_list[index] ⊇ tmp_den_collection ||
          check_Pak_covering(
            den_collection_list[index], tmp_den_collection;
            SP_replacements=kinematic_relations
          ),
        remaining_indices, index_index + 1
      )
      isnothing(tmp_index_index) && break
      index_index, largest_index = tmp_index_index, remaining_indices[tmp_index_index]
    end

    key_collection = den_collection_list[largest_index]
    repl_rules_list = Vector{Union{Dict{Basic, Basic}, Nothing}}(undef, length(remaining_indices))
    for (ii, index) ∈ enumerate(remaining_indices)
      target_collection = den_collection_list[index]
      if index == largest_index || key_collection ⊇ target_collection
        repl_rules_list[ii] = Dict{Basic, Basic}() # do not check the same topology or directly sub-topology
        continue
      end # if
      @info "Checking $largest_index covering: $(ii) / $(length(remaining_indices))."
      repl_rules_list[ii] = find_Pak_momentum_shifts(
        key_collection, target_collection;
        SP_replacements=kinematic_relations,
        find_external_momentum_shifts=find_external_momentum_shifts
      ) # end find_Pak_momentum_shifts
    end # for index ∈ remaining_indices
      
    indices = findall(!isnothing, repl_rules_list)
    push!(result_den_collection_list, key_collection)
    push!(result_repl_rules_list,
      Dict{Int, Dict{Basic, Basic}}(remaining_indices[ii] => repl_rules_list[ii] for ii ∈ indices)
    ) # end push!
    deleteat!(remaining_indices, indices)
    # print("\r")
  end # while

  @info "Found $(length(result_den_collection_list)) topologies."

  return result_den_collection_list, result_repl_rules_list
end # function construct_den_topology
