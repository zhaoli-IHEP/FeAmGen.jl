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

  independent_external_momenta = reduce(union, map(dc -> dc.external_momenta, den_collection_list))
  # independent_external_momenta = external_momenta[begin:end-1]
  independent_external_momenta = map(Basic, independent_external_momenta)
  n_ind_ext = length(independent_external_momenta)

  n_sp = binomial(n_loop, 2) + n_loop * (n_ind_ext + 1)

  den_list_list = map(dc -> unique(dc.denominators), den_collection_list)

  graded_indices_list = [ [ [ii] for ii ∈ eachindex(den_list_list) ] ]
  previous_indices_list = last(graded_indices_list)

  n_den_collect = length(den_list_list)
  counter = n_den_collect - 1
  while counter > 0
    this_indices_list = Vector{Int}[]
    for one_indices ∈ previous_indices_list
      for one_index ∈ (last(one_indices) + 1):n_den_collect
        tmp_indices = vcat(one_indices, [one_index])
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

  @assert (isempty ∘ symdiff)(reduce(union, final_indices_list), eachindex(den_list_list))

  result_den_collection_list = FeynmanDenominatorCollection[]
  for indices ∈ final_indices_list
    den_list = reduce(union, den_list_list[indices])
    push!(result_den_collection_list,
      FeynmanDenominatorCollection(
        loop_momenta, independent_external_momenta, den_list;
        check_validity=false
      )
    ) # end push!
  end # for indices

  return result_den_collection_list, final_indices_list
end # function minimize_topology_list_directly

function make_complete_topology(
  ::Val{:Canonicalization},
  topology::FeynmanDenominatorCollection
)::FeynmanDenominatorCollection
  n_loop = loop_momenta_number(topology)
  @assert n_loop ∈ 1:4 "Only support 1-4 loop(s), but got $(n_loop) loop(s)!"
  loop_momenta = map(Basic, topology.loop_momenta)
  q_list = [Basic("q$ii") for ii ∈ 1:n_loop]
  external_momenta = topology.external_momenta
  k_list = map(Basic, external_momenta#=[begin:end-1]=#)
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
end # function make_complete_topology

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
  k_list = map(Basic, denominator_collection.external_momenta#=[begin:end-1]=#)
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
