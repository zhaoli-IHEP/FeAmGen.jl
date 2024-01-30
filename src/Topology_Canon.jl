function minimize_topology_list_directly(
  ::Val{:Canonicalization},
  den_collection_list::Vector{FeynmanDenominatorCollection}
)::Vector{FeynmanDenominatorCollection}
  # initial information #######################################################
  loop_momenta = reduce(union, map(dc -> dc.loop_momenta, den_collection_list))
  n_loop = length(loop_momenta)

  preferred_loop_exprs_list = if haskey(preferred_vac_mom_dict(), n_loop)
    exprs_list = preferred_vac_mom_dict()[n_loop]
    for (ii, exprs) ∈ enumerate(exprs_list)
      map!(x -> expand(x^2), exprs_list[ii], exprs)
    end # for
    exprs_list
  else
    @assert false "No preferred loop expressions for $n_loop-loop! Please use `:PakAlgorithm` instead!"
  end # if

  independent_external_momenta = reduce(union, map(dc -> dc.external_momenta, den_collection_list))
  pop!(independent_external_momenta)
  n_ind_ext = length(independent_external_momenta)

  n_sp = binomial(n_loop, 2) + n_loop * n_ind_ext

  den_list_list = map(dc -> dc.denominators, den_collection_list)
  # end initial information ###################################################

  graded_indices_list = [ [ [ii] for ii ∈ eachindex(den_collection_list) ] ]
  previous_indices_list = last(graded_indices_list)
  previous_union_den_list = map(
    indices -> den_collection_list[indices],
    previous_indices_list
  )

  n_den_collect = length(den_collection_list)
  counter = n_den_collect - 1
  while counter > 0
    this_indices_list = Vector{Int}[]
    this_union_den_list = Vector{FeynmanDenominator}[]
    for (ii, one_indices) ∈ enumerate(previous_indices_list)
      for one_index ∈ (last(one_indices) + 1):n_den_collect
        push!(this_indices_list, (sort! ∘ union)(one_indices, [one_index]))
        push!(this_union_den_list,
          union(previous_union_den_list[ii], den_list_list[[one_index]])
        )
      end # for one_index
    end # for (ii, one_indices)
    
    valid_indices_positions = findall(
      union_den -> length(union_den) ≤ n_sp &&
        calculate_denominator_collection_rank(
          loop_momenta, independent_external_momenta, union_den
        ) == length(union_den),
      this_union_den_list
    ) # end valid_indices_positions
    isempty(valid_indices_positions) && break
    this_indices_list = this_indices_list[valid_indices_positions]
    this_union_den_list = this_union_den_list[valid_indices_positions]

    valid_indices_positions = Int[]
    for (ii, union_den) ∈ enumerate(this_union_den_list)
      pure_loop_exprs = map(get_denominator_expr, union_den)
      map!(expr -> expand(expr - subs(expr, Dict(loop_momenta .=> 0))),
        pure_loop_exprs, pure_loop_exprs
      ) # end map!
      unique!(pure_loop_exprs)
      any(preferred_loop_exprs -> pure_loop_exprs ⊆ preferred_loop_exprs,
        preferred_loop_exprs_list
      ) || continue
      push!(valid_indices_positions, ii)
    end # for (ii, union_den)
    isempty(valid_indices_positions) && break
    this_indices_list = this_indices_list[valid_indices_positions]
    this_union_den_list = this_union_den_list[valid_indices_positions]

    push!(graded_indices_list, this_indices_list)
    previous_indices_list = this_indices_list
    previous_union_den_list = this_union_den_list

    counter -= 1
  end # while

  reverse!(graded_indices_list)
  final_indices_list = Vector{Int}[]
  for indices_list ∈ graded_indices_list
    covered_indices = reduce(union, final_indices_list)
    to_be_added_indices_list = Vector{Int}[]
    for indices ∈ indices_list
      !isempty(indices ∩ covered_indices) && continue
      push!(to_be_added_indices_list, indices)
    end # for indices
    to_be_added_indices_list = get_greedy_cover(to_be_added_indices_list)
    for indices ∈ to_be_added_indices_list
      push!(final_indices_list, indices)
    end # for indices
  end # for indices_list

  result_den_collection_list = FeynmanDenominatorCollection[]
  for indices ∈ final_indices_list
    push!(
      result_den_collection_list,
      reduce(union, den_collection_list[indices])
    ) # end push!
  end # for indices

  return result_den_collection_list
end # function minimize_topology_list_directly
