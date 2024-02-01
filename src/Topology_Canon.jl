function minimize_topology_list_directly(
  ::Val{:Canonicalization},
  den_collection_list::Vector{FeynmanDenominatorCollection}
)::Tuple{
  Vector{FeynmanDenominatorCollection}, # minimal topologies
  Vector{Vector{Int}} # covering indices
}
  # initial information #######################################################
  loop_momenta = reduce(union, map(dc -> dc.loop_momenta, den_collection_list))
  loop_momenta = map(Basic, loop_momenta)
  n_loop = length(loop_momenta)

  preferred_loop_exprs_list = if haskey(preferred_vac_mom_dict(), n_loop)
    exprs_list = preferred_vac_mom_dict()[n_loop]
    for (ii, exprs) ∈ enumerate(exprs_list)
      map!(x -> expand(x^2), exprs_list[ii], exprs)
    end # for
    exprs_list
  else
    error("No preferred loop expressions for $n_loop-loop! Please use `:PakAlgorithm` instead!")
  end # if

  external_momenta = reduce(union, map(dc -> dc.external_momenta, den_collection_list))
  independent_external_momenta = external_momenta[begin:end-1]
  independent_external_momenta = map(Basic, independent_external_momenta)
  n_ind_ext = length(independent_external_momenta)

  n_sp = binomial(n_loop, 2) + n_loop * (n_ind_ext + 1)

  original_den_list_list = map(dc -> unique(dc.denominators), den_collection_list)
  # end initial information ###################################################

  den_list_list = Vector{FeynmanDenominator}[]
  original_den_index_list = (collect ∘ eachindex)(original_den_list_list)
  covering_indices = Vector{Int}[]
  while !isempty(original_den_list_list)
    _, index = findmax(length, original_den_list_list)
    target_den_list = original_den_list_list[index]
    push!(den_list_list, target_den_list)

    tmp_indices = findall(den_list -> den_list ⊆ target_den_list, original_den_list_list)
    push!(covering_indices, original_den_index_list[tmp_indices])
    deleteat!(original_den_list_list, tmp_indices)

    filter!(den_list -> den_list ⊈ target_den_list, original_den_list_list)
  end # while

  graded_indices_list = [ [ [ii] for ii ∈ eachindex(den_list_list) ] ]
  previous_indices_list = last(graded_indices_list)

  n_den_collect = length(den_list_list)
  counter = n_den_collect - 1
  while counter > 0
    # @show counter
    this_indices_list = Vector{Int}[]
    # this_union_den_list = Vector{FeynmanDenominator}[]
    for one_indices ∈ previous_indices_list
      for one_index ∈ (last(one_indices) + 1):n_den_collect
        tmp_indices = (sort! ∘ union)(one_indices, [one_index])
        tmp_union_den = reduce(union, den_list_list[tmp_indices])
        length(tmp_union_den) ≤ n_sp || continue
        calculate_denominator_collection_rank(
          loop_momenta, independent_external_momenta, tmp_union_den
        ) == length(tmp_union_den) || continue

        # :Canonicalization specific ##########################################
        pure_loop_exprs = Vector{Basic}(undef, length(tmp_union_den))
        for (ii, den) ∈ enumerate(tmp_union_den)
          momentum_set = den.denominator_momenta
          @assert length(momentum_set) == 1 "Not a standard Feynman propagator denominator!"
          pure_loop_exprs[ii] = (Basic ∘ first)(momentum_set)
        end # for (ii, den)
        map!(expr -> expand(expr - subs(expr, Dict(loop_momenta .=> zero(Basic)))),
          pure_loop_exprs, pure_loop_exprs
        ) # end map!
        unique!(pure_loop_exprs)
        length(pure_loop_exprs) ≤ 3 * (n_loop - 1) || continue
        any(preferred_loop_exprs -> pure_loop_exprs ⊆ preferred_loop_exprs,
          preferred_loop_exprs_list
        ) || continue
        # end :Canonicalization specific ######################################

        push!(this_indices_list, tmp_indices)
      end # for one_index
    end # for (ii, one_indices)

    isempty(this_indices_list) && break
    push!(graded_indices_list, this_indices_list)
    previous_indices_list = this_indices_list

    counter -= 1
  end # while

  final_indices_list = (get_greedy_cover ∘ pop!)(graded_indices_list)
  while !isempty(graded_indices_list)
    # @show final_indices_list
    tmp_indices_list = pop!(graded_indices_list)
    covered_indices = reduce(union, final_indices_list)
    to_be_added_indices_list = Vector{Int}[]
    for indices ∈ tmp_indices_list
      isempty(indices ∩ covered_indices) || continue
      push!(to_be_added_indices_list, indices)
    end # for indices
    isempty(to_be_added_indices_list) && continue
    to_be_added_indices_list = get_greedy_cover(to_be_added_indices_list)
    final_indices_list = append!(final_indices_list, to_be_added_indices_list)
  end

  # @show final_indices_list
  @assert (isempty ∘ symdiff)(reduce(union, final_indices_list), eachindex(den_list_list))

  result_den_collection_list = FeynmanDenominatorCollection[]
  result_covering_indices = Vector{Int}[]
  for indices ∈ final_indices_list
    den_list = reduce(union, den_list_list[indices])
    push!(result_den_collection_list,
      FeynmanDenominatorCollection(
        loop_momenta, external_momenta, den_list;
        check_validity=false
      )
    ) # end push!
    push!(result_covering_indices, reduce(union, covering_indices[indices]))
  end # for indices

  return result_den_collection_list, result_covering_indices
end # function minimize_topology_list_directly

function construct_den_topology(
  ::Val{:Canonicalization},
  den_collection_list::Vector{FeynmanDenominatorCollection}
)::Tuple{
  Vector{FeynmanDenominatorCollection},
  Vector{Dict{Int, Dict{Basic, Basic}}}
}
  result_den_collection_list = FeynmanDenominatorCollection[]
  result_repl_rules_list = Dict{Int, Dict{Basic, Basic}}[]
  remaining_indices = (collect ∘ eachindex)(den_collection_list)

  @info "Minimalizing the topologies with Canonicalization algorithm."
  while !isempty(remaining_indices)
    @info "Number of remaining topologies: $(length(remaining_indices))."

    _, index_index = findmax(length, den_collection_list[remaining_indices])
    largest_index = remaining_indices[index_index]
    while true
      tmp_den_collection = den_collection_list[largest_index]
      tmp_index_index = findnext(
        index -> tmp_den_collection ≠ den_collection_list[index] &&
          !is_Canonicalization_equivalent(
            tmp_den_collection, den_collection_list[index]
          ) && den_collection_list[index] ⊇ tmp_den_collection ||
          check_Canonicalization_covering(
            den_collection_list[index], tmp_den_collection
          ),
        remaining_indices, index_index + 1
      )
      isnothing(tmp_index_index) && break
      index_index, largest_index = tmp_index_index, remaining_indices[tmp_index_index]
    end # while

    key_collection = den_collection_list[largest_index]
    repl_rules_list = Vector{Union{Dict{Basic, Basic}, Nothing}}(undef, length(remaining_indices))
    for (ii, index) ∈ enumerate(remaining_indices)
      target_collection = den_collection_list[index]
      if index == largest_index || key_collection ⊇ target_collection
        repl_rules_list[ii] = Dict{Basic, Basic}() # do not check the same topology or directly sub-topology
        continue
      end # if
      @info "Checking $largest_index covering: $(ii) / $(length(remaining_indices))."
      repl_rules_list[ii] = find_Canonicalization_momentum_shifts(
        key_collection, target_collection
      ) # end find_Canonicalization_momentum_shifts
    end # for (ii, index)

    indices = findall(!isnothing, repl_rules_list)
    push!(result_den_collection_list, key_collection)
    push!(result_repl_rules_list,
      Dict{Int, Dict{Basic, Basic}}(remaining_indices[ii] => repl_rules_list[ii] for ii ∈ indices)
    ) # end push!
    deleteat!(remaining_indices, indices)
  end # while

  @info "Constructed $(length(result_den_collection_list)) topologies."

  return result_den_collection_list, result_repl_rules_list
end # function construct_den_topology

function make_complete_topology(
  ::Val{:Canonicalization},
  topology::FeynmanDenominatorCollection
)::FeynmanDenominatorCollection
  n_loop = loop_momenta_number(topology)
  @assert n_loop ∈ 1:4 "Only support 1-4 loop(s), but got $(n_loop) loop(s)!"
  loop_momenta = map(Basic, topology.loop_momenta)
  q_list = [Basic("q$ii") for ii ∈ 1:n_loop]
  external_momenta = topology.external_momenta
  k_list = map(Basic, external_momenta[begin:end-1])
  n_ind_ext = length(k_list)
  tmp_dict = Dict(loop_momenta .=> q_list)

  n_sp = binomial(n_loop, 2) + n_loop * (n_ind_ext + 1)

  denominator_momenta = Vector{Basic}(undef, length(topology))
  for (ii, den) ∈ enumerate(topology)
    momentum_set = den.denominator_momenta
    @assert length(momentum_set) == 1 "Not a standard Feynman propagator denominator!"
    denominator_momenta[ii] = (Basic ∘ first)(momentum_set)
  end # for (ii, den)
  map!(den_mom -> subs(den_mom, tmp_dict), denominator_momenta, denominator_momenta)

  new_denominator_list = Vector{FeynmanDenominator}(undef, length(topology))
  for (ii, den) ∈ enumerate(topology)
    new_denominator_list[ii] = FeynmanDenominator(
      denominator_momenta[ii], den.mass, den.width;
      loop_momenta=q_list,
      external_momenta=den.external_momenta,
      check_validity=false
    ) # end FeynmanDenominator
  end # for (ii, den)

  vacuum_denominator_momenta = begin
    preferred_loop_exprs_list = preferred_vac_mom_dict()[n_loop]
    pure_loop_exprs = map(expr -> expand(expr - subs(expr, Dict(q_list .=> 0))),
      denominator_momenta
    ) # end map
    unique!(pure_loop_exprs)
    selected_index = findfirst(
      preferred_loop_exprs -> pure_loop_exprs ⊆ preferred_loop_exprs,
      preferred_loop_exprs_list
    ) # end findfirst
    @assert !isnothing(selected_index)
    preferred_loop_exprs_list[selected_index]
  end

  for k ∈ vcat(zeros(Basic, 1), k_list), p ∈ vacuum_denominator_momenta, relative_sign ∈ [1, -1]
    trial_denominator = FeynmanDenominator(
      p + relative_sign * k, 0, 0;
      loop_momenta=q_list,
      external_momenta=external_momenta,
      check_validity=false
    )
    trial_denominator ∈ topology && continue
    trial_denominator_list = union(new_denominator_list, [trial_denominator])
    the_rank = calculate_denominator_collection_rank(
      loop_momenta, k_list, trial_denominator_list
    )
    the_rank == length(trial_denominator_list) || continue
    new_denominator_list = trial_denominator_list
    the_rank == n_sp && break
  end # for k, p, relative_sign

  return FeynmanDenominatorCollection(
    q_list, external_momenta, new_denominator_list;
    check_validity=false
  ) # end FeynmanDenominatorCollection
end

function is_Canonicalization_equivalent(
  d1::FeynmanDenominatorCollection,
  d2::FeynmanDenominatorCollection
)::Bool
  length(d1) == length(d2) || return false
  loop_momenta_number(d1) == loop_momenta_number(d2) || return false

  canonical_d1 = canonicalize_denominator_collection(d1)
  canonical_d2 = canonicalize_denominator_collection(d2)

  return canonical_d1 == canonical_d2
end # function is_Canonicalization_equivalent

check_Canonicalization_covering(
  d1::FeynmanDenominatorCollection,
  d2::FeynmanDenominatorCollection
)::Bool = (!isnothing ∘ find_Canonicalization_momentum_shifts)(d1, d2)

function find_Canonicalization_momentum_shifts(
  d1::FeynmanDenominatorCollection,
  d2::FeynmanDenominatorCollection
)::Union{Dict{Basic, Basic}, Nothing}
  n_d1 = length(d1)
  n_d2 = length(d2)
  n_d1 ≥ n_d2 || return nothing
  if loop_momenta_number(d1) ≠ loop_momenta_number(d2)
    @warn "Do not support different loop momenta number: ($(loop_momenta_number(d1_collection)) vs $(loop_momenta_number(d2_collection)))."
    return nothing
  end # if
  if external_momenta_number(d1) ≠ external_momenta_number(d2)
    @warn "Do not support different external momenta number: ($(external_momenta_number(d1_collection)) vs $(external_momenta_number(d2_collection)))."
    return nothing
  end # if

  canonical_d2, d2_momentum_shifts = canonicalize_denominator_collection(d2; return_momentum_shifts=true)
  for selected_indices ∈ combinations(1:n_d1, n_d2)
    sub_d1 = get_sub_collection(d1, selected_indices)
    canonical_sub_d1 = canonicalize_denominator_collection(sub_d1)
    canonical_sub_d1 == canonical_d2 && return d2_momentum_shifts
  end # for selected_indices

  return nothing
end # function find_Canonicalization_momentum_shifts

function canonicalize_denominator_collection(
  denominator_collection::FeynmanDenominatorCollection;
  return_momentum_shifts::Bool=false
)::Union{
  FeynmanDenominatorCollection,
  Tuple{
    FeynmanDenominatorCollection,
    Dict{Basic, Basic}
  }
}
  n_loop = loop_momenta_number(denominator_collection)
  @assert n_loop ∈ 1:4 "Only support 1-4 loop(s), but got $(n_loop) loop(s)!"
  loop_momenta = map(Basic, denominator_collection.loop_momenta)
  q_list = [Basic("q$ii") for ii ∈ 1:n_loop]
  k_list = map(Basic, denominator_collection.external_momenta[begin:end-1])
  tmp_dict = Dict(loop_momenta .=> q_list)

  denominator_momenta = Vector{Basic}(undef, length(denominator_collection))
  for (ii, den) ∈ enumerate(denominator_collection)
    momentum_set = den.denominator_momenta
    @assert length(momentum_set) == 1 "Not a standard Feynman propagator denominator!"
    denominator_momenta[ii] = (Basic ∘ first)(momentum_set)
  end # for (ii, den)
  map!(den_mom -> subs(den_mom, tmp_dict), denominator_momenta, denominator_momenta)

  canonicalization_map = gen_loop_mom_canon_map(denominator_momenta, q_list, k_list)
  map!(den_mom -> subs(den_mom, canonicalization_map), denominator_momenta, denominator_momenta)

  new_denominator_list = Vector{FeynmanDenominator}(undef, length(denominator_collection))
  for (ii, den) ∈ enumerate(denominator_collection)
    new_denominator_list[ii] = FeynmanDenominator(
      denominator_momenta[ii], den.mass, den.width;
      loop_momenta=q_list,
      external_momenta=den.external_momenta,
      check_validity=false
    ) # end FeynmanDenominator
  end # for (ii, den)

  canonical_denominator_collection = FeynmanDenominatorCollection(
    q_list,
    denominator_collection.external_momenta,
    new_denominator_list;
    check_validity=false
  ) # end FeynmanDenominatorCollection
  return return_momentum_shifts ? (
    canonical_denominator_collection,
    canonicalization_map
  ) : canonical_denominator_collection
end # function canonicalize_denominator_collection
