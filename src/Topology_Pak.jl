function construct_topology(
  # ::Val{:PakAlgorithm},
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
    largest_den_collection = den_collection_list[largest_index]
    # for ii ∈ (index_index+1):length(remaining_indices)
    #   @info "Checking the largest index: $ii / $(length(remaining_indices))."
    #   index = remaining_indices[ii]
    #   tmp_den_collection = den_collection_list[index]
    #   tmp_den_collection == largest_den_collection && continue
    #   is_Pak_equivalent(
    #     tmp_den_collection, largest_den_collection;
    #     SP_replacements=kinematic_relations
    #   ) && continue
    #   tmp_den_collection ⊇ largest_den_collection ||
    #     check_Pak_covering(
    #       tmp_den_collection, largest_den_collection;
    #       SP_replacements=kinematic_relations
    #     ) || continue
    #   @info "Found a larger index: $ii / $(length(remaining_indices))."
    #   largest_index = index
    #   largest_den_collection = tmp_den_collection
    # end # for ii

    repl_rules_list = Vector{Union{Dict{Basic, Basic}, Nothing}}(undef, length(remaining_indices))
    for (ii, index) ∈ enumerate(remaining_indices)
      target_collection = den_collection_list[index]
      if index == largest_index || largest_den_collection ⊇ target_collection
        repl_rules_list[ii] = Dict{Basic, Basic}() # do not check the same topology or directly sub-topology
        continue
      end # if
      @info "Checking $largest_index covering: $(ii) / $(length(remaining_indices))."
      repl_rules_list[ii] = find_Pak_momentum_shifts(
        largest_den_collection, target_collection;
        SP_replacements=kinematic_relations,
        find_external_momentum_shifts=find_external_momentum_shifts
      ) # end find_Pak_momentum_shifts
    end # for (ii, index)
      
    positions = findall(!isnothing, repl_rules_list)
    push!(result_den_collection_list, largest_den_collection)
    push!(result_repl_rules_list,
      Dict{Int, Dict{Basic, Basic}}(remaining_indices[ii] => repl_rules_list[ii] for ii ∈ positions)
    ) # end push!
    deleteat!(remaining_indices, positions)
  end # while

  @info "Constructed $(length(result_den_collection_list)) topologies."

  return result_den_collection_list, result_repl_rules_list
end # function construct_den_topology
