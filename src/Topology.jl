include("Topology_legacy.jl")
include("Topology_Canon.jl")
include("Topology_Pak.jl")

function construct_den_topology(
  amp_dir::String;
  mode::Symbol=:Canonicalization, # :Canonicalization or :PakAlgorithm
  check_momentum_shifts::Bool=true,
  # use_reference_topology::Bool=false,
  # reference_topology_directory::String="",
  # reference_topologies::Vector{FeynmanDenominatorCollection}=FeynmanDenominatorCollection[],
  return_complete_topology::Bool=true,
  check_all_kinematic_relations::Bool=false,
  find_external_momentum_shifts::Bool=false,
  recheck_momentum_shifts::Bool=false
)::Vector{FeynmanDenominatorCollection}

  # check options #############################################################
  @assert mode ∈ [
    :Canonicalization,
    :PakAlgorithm
  ] "Do not support mode $mode."

  # use_reference_topology &&
  #   return construct_den_topology_with_reference_topology(
  #     amp_dir;
  #     mode=mode,
  #     check_momentum_shifts=check_momentum_shifts,
  #     reference_topology_directory=reference_topology_directory,
  #     reference_topologies=reference_topologies,
  #     check_all_kin_relation=check_all_kin_relation,
  #     find_external_momentum_shifts=find_external_momentum_shifts
  #   )
  # end check options #########################################################

  # init ######################################################################
  amp_filename_list, amp_den_collect_list, kin_relation_dict  = read_loop_denominators(
    Val(:AmplitudeDirectory), amp_dir, check_all_kinematic_relations
  )
  covering_dict = Dict{FeynmanDenominatorCollection, Dict{Int, Dict{Basic, Basic}}}()

  top_dir_list = splitpath(amp_dir)
  topology_name = if mode == :Canonicalization
    "Canonicalization_topologies"
  elseif mode == :PakAlgorithm
    "Pak_topologies"
  end # if
  top_dir_list[end] = if endswith(top_dir_list[end], "amplitudes")
    replace(top_dir_list[end], "amplitudes" => topology_name)
  else
    top_dir_list[end] *= "_" * topology_name
  end
  topology_directory = joinpath(top_dir_list)
  bk_mkdir(topology_directory)
  # end init ##################################################################

  # step 0 ####################################################################
  # ⊆, == #####################################################################
  @info "Removing redundant topologies."
  init_den_collect_list = FeynmanDenominatorCollection[]
  tmp_den_collect_list = copy(amp_den_collect_list)
  tmp_den_collect_indices = (collect ∘ eachindex)(tmp_den_collect_list)
  init_covering_indices = Vector{Int}[]
  while !isempty(tmp_den_collect_list)
    _, index = findmax(length, tmp_den_collect_list)
    to_be_chosen_den_collect = tmp_den_collect_list[index]

    to_be_chosen_den_collect_list = filter(
      den_collect -> length(den_collect) == length(to_be_chosen_den_collect),
      tmp_den_collect_list
    ) # end filter
    while true
      neq_to_be_chosen_den_collect_list = filter(
        den_collect -> den_collect ≠ to_be_chosen_den_collect,
        to_be_chosen_den_collect_list
      ) # end filter
      filter!(
        den_collect -> den_collect ⊇ to_be_chosen_den_collect,
        neq_to_be_chosen_den_collect_list
      ) # end filter!
      isempty(neq_to_be_chosen_den_collect_list) && break
      to_be_chosen_den_collect = popfirst!(neq_to_be_chosen_den_collect_list)
    end # while

    push!(init_den_collect_list, to_be_chosen_den_collect)

    to_be_removed_positions = findall(
      den_collect -> den_collect ⊆ to_be_chosen_den_collect,
      tmp_den_collect_list
    ) # end findall
    push!(init_covering_indices, tmp_den_collect_indices[to_be_removed_positions])
    deleteat!(tmp_den_collect_indices, to_be_removed_positions)
    deleteat!(tmp_den_collect_list, to_be_removed_positions)
  end # while
  @info "Got $(length(init_den_collect_list)) initial topologies."
  # end step 0 ################################################################

  if !check_momentum_shifts
    @info "Minimizing the initial topologies without checking momentum shifts."
    greedy_den_collect_list, greedy_covering_indices = minimize_topology_list_directly(
      Val(mode), init_den_collect_list
    ) # end minimize_topology_list_directly
    @info "Done and got $(length(greedy_den_collect_list)) topologies."

    complete_topologies = greedy_den_collect_list
    if return_complete_topology
      @info "Constructing complete topologies."
      complete_topologies = if mode == :Canonicalization
        make_complete_topology.(Val(:Canonicalization), greedy_den_collect_list)
      else
        @warn "Do not support to complete topologies with Pak algorithm yet."
        greedy_den_collect_list
      end # if
    end # if

    for (ii, topology) ∈ enumerate(complete_topologies)
      all_covering_indices = reduce(union, init_covering_indices[greedy_covering_indices[ii]])
      covering_dict[topology] = Dict{Int, Dict{Basic, Basic}}(all_covering_indices .=> Ref(Dict{Basic, Basic}()))
    end # for (ii, topology)

    write_topology(covering_dict, amp_filename_list, topology_directory, kin_relation_dict)

    return complete_topologies
  end # if

  # step 1 ####################################################################
  # find momentum shifts #######################################################
  first_shifted_den_collect_list, first_shifted_repl_rules_list =
    construct_topology(
      init_den_collect_list,
      kin_relation_dict,
      find_external_momentum_shifts
    ) # end construct_den_topology
  # end step 1 ################################################################

  # step 2 ####################################################################
  # greedy minimize topologies ################################################
  @info "Minimalizing the initial topologies."
  greedy_den_collect_list, greedy_covering_indices = minimize_topology_list_directly(
    Val(mode), first_shifted_den_collect_list
  ) # end minimal_topology_list
  @info "Done and got $(length(greedy_den_collect_list)) topologies."
  # end step 2 ################################################################

  # step 3 ####################################################################
  # find momentum shifts again ################################################
  second_shifted_den_collect_list, second_shifted_repl_rules_list =
    construct_topology(
      greedy_den_collect_list,
      kin_relation_dict,
      find_external_momentum_shifts
    ) # end construct_den_topology
  # end step 3 ################################################################

  # step 4 ####################################################################
  # complete topologies #######################################################
  complete_topologies = second_shifted_den_collect_list
  if return_complete_topology
    @info "Constructing complete topologies."
    complete_topologies = if mode == :Canonicalization
      make_complete_topology.(Val(:Canonicalization), second_shifted_den_collect_list)
    else
      @warn "Do not support to complete topologies with (deep) Pak algorithm yet."
      second_shifted_den_collect_list
    end # if
  end # if
  # end step 4 ################################################################

  # construct covering_dict ###################################################
  covering_dict = Dict{FeynmanDenominatorCollection, Dict{Int, Dict{Basic, Basic}}}()
  if recheck_momentum_shifts
    @info "Re-checking all original topologies are covered by the constructed topologies."
    for (ii, topology) ∈ enumerate(complete_topologies)
      repl_rules = Dict{Int, Dict{Basic, Basic}}()
      for (jj, den_collection) ∈ enumerate(amp_den_collect_list)
        @info "Checking if the constructed topology #$ii (total: $(length(result_topology_list))) covering the original topology #$jj (total: $(length(amp_den_collect_list)))."
        if mode == :Canonicalization
          momentum_shifts = find_Canonicalization_momentum_shifts(
            topology, den_collection
          ) # end find_Canonicalization_momentum_shifts
        elseif mode == :PakAlgorithm
          momentum_shifts = find_Pak_momentum_shifts(
            topology, den_collection;
            SP_replacements=kin_relation_dict,
            find_external_momentum_shifts=find_external_momentum_shifts
          ) # end find_Pak_momentum_shifts
          isnothing(momentum_shifts) && continue
          repl_rules[jj] = momentum_shifts
        end # if
      end # for (jj, den_collection)
      covering_dict[topology] = repl_rules
    end # for (ii, topology)
  else
    for (ii, topology) ∈ enumerate(complete_topologies)
      tmp_dict = Dict{Int, Dict{Basic, Basic}}()
      for (second_key, second_momentum_shifts) ∈ second_shifted_repl_rules_list[ii]
        for greedy_covering_index ∈ greedy_covering_indices[second_key]
          for (first_key, first_momentum_shifts) ∈ first_shifted_repl_rules_list[greedy_covering_index]
            final_momentum_shifts = if isempty(first_momentum_shifts)
              second_momentum_shifts
            elseif isempty(second_momentum_shifts)
              first_momentum_shifts
            else
              tmp_momentum_shifts = Dict{Basic, Basic}()
              for (key, value) ∈ first_momentum_shifts
                tmp_momentum_shifts[key] = (expand ∘ subs)(value, second_momentum_shifts)
              end # for (key, value)
              tmp_momentum_shifts
            end # if

            for amp_index ∈ init_covering_indices[first_key]
              tmp_dict[amp_index] = final_momentum_shifts
            end # for amp_index
          end # for (first_key, first_momentum_shifts)
        end # for greedy_covering_index
      end # for (second_key, second_momentum_shifts)
      covering_dict[topology] = tmp_dict
    end # for (ii, topology)
  end # if

  @assert (isempty ∘ symdiff)(
    reduce(union, map(collect ∘ keys, (collect ∘ values)(covering_dict))),
    eachindex(amp_den_collect_list)
  )
  # end construct covering_dict ###############################################

  # write topology ############################################################
  write_topology(covering_dict, amp_filename_list, topology_directory, kin_relation_dict)
  # end write topology ########################################################

  return complete_topologies
end # function construct_topology

function read_loop_denominators(
  ::Val{:AmplitudeDirectory},
  amplitude_directory::String,
  check_all_kinematic_relations::Bool
)::Tuple{
  Vector{String}, # amplitude file list
  Vector{FeynmanDenominatorCollection},
  Dict{Basic, Basic} # kinematic relations
}
  @assert isdir(amplitude_directory) "Got non-existent directory: $amplitude_directory."
  amp_filename_list = filter!(endswith(".jld2"), readdir(amplitude_directory; join=true, sort=false))
  sort!(amp_filename_list, by=filename -> begin
                                            main_name, _ = (splitext ∘ basename)(filename)
                                            index_str = match(r"[1-9]\d*", main_name).match
                                            parse(Int, index_str)
                                          end)
  @assert !isempty(amp_filename_list) "Cannot find any amplitude file with extension .jld2 in $amplitude_directory."
  first_den_collection, kin_relation_dict = read_loop_denominators(
    Val(:AmplitudeFile), first(amp_filename_list), true
  )
  den_collections = FeynmanDenominatorCollection[first_den_collection]
  for amp_file ∈ amp_filename_list[2:end]
    den_collection, new_kin_relation_dict = read_loop_denominators(
      Val(:AmplitudeFile), amp_file, check_all_kinematic_relations
    )
    push!(den_collections, den_collection)
    isempty(new_kin_relation_dict) && continue
    for (key, value) ∈ new_kin_relation_dict
      if haskey(kin_relation_dict, key)
        @assert kin_relation_dict[key] == value """
        \rGot different kinematic relation:
        \r    - $key => $(kin_relation_dict[key]) in previous files;
        \r    - $key => $value in $amp_file.
        """
      else
        kin_relation_dict[key] = value
      end # if
    end # for (key, value)
  end # for amp_file
  return amp_filename_list, den_collections, kin_relation_dict
end # function read_loop_denominators
function read_loop_denominators(
  ::Val{:AmplitudeFile},
  amplitude_file_name::String,
  read_kinematic_relations::Bool
)::Tuple{
  FeynmanDenominatorCollection,
  Dict{Basic, Basic} # kinematic relations
}
  amplitude_file = load(amplitude_file_name)
  n_loop = amplitude_file["n_loop"]
  loop_den_list = map(Basic, amplitude_file["loop_den_list"])
  loop_momenta = ["q$ii" for ii ∈ 1:n_loop]
  external_momenta = subs.(
    map(Basic, amplitude_file["ext_mom_list"]),
    (Ref ∘ Dict{Basic, Basic})(Basic(key) => Basic(value) for (key, value) ∈ amplitude_file["mom_symmetry"])
  ) # end external_momenta
  unique!(external_momenta)
  den_list = Vector{FeynmanDenominator}(undef, length(loop_den_list))
  for (den_index, loop_den) ∈ enumerate(loop_den_list)
    @assert get_name(loop_den) == "Den"
    mom, mass, width = get_args(loop_den)
    den_list[den_index] = FeynmanDenominator(mom, mass, width)
  end # for

  kin_relation_dict = Dict{Basic, Basic}()
  read_kinematic_relations ||
    return FeynmanDenominatorCollection(loop_momenta, external_momenta, den_list), kin_relation_dict

  for (key, value) ∈ amplitude_file["kin_relation"]
    key = Basic(key)
    value = Basic(value)
    get_name(key) == "SP" || continue
    kin_relation_dict[key] = value
  end # for (key, value)

  return FeynmanDenominatorCollection(loop_momenta, external_momenta, den_list), kin_relation_dict
end # function read_loop_denominators

function minimize_topology_list_directly(
  ::Any,
  den_collection_list::Vector{FeynmanDenominatorCollection}
)::Tuple{
  Vector{FeynmanDenominatorCollection}, # minimal topologies
  Vector{Vector{Int}} # covering indices
}
  # initial information #######################################################
  loop_momenta = reduce(union, map(dc -> dc.loop_momenta, den_collection_list))
  loop_momenta = map(Basic, loop_momenta)
  n_loop = length(loop_momenta)

  external_momenta = reduce(union, map(dc -> dc.external_momenta, den_collection_list))
  independent_external_momenta = external_momenta[begin:end-1]
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
        loop_momenta, external_momenta, den_list;
        check_validity=false
      )
    ) # end push!
  end # for indices

  return result_den_collection_list, final_indices_list
end # function minimize_topology_list_directly

function is_complete_topology(den_collection::FeynmanDenominatorCollection)::Bool
  loop_momenta = map(Basic, den_collection.loop_momenta)
  n_loop = length(loop_momenta)
  independent_external_momenta = map(Basic, den_collection.external_momenta[begin:end-1])
  n_ind_ext = length(independent_external_momenta)
  n_sp = binomial(n_loop, 2) + n_loop * (n_ind_ext + 1)

  return calculate_denominator_collection_rank(
    loop_momenta, independent_external_momenta, den_collection.denominators
  ) == n_sp
end # function is_complete_topology

function calculate_denominator_collection_rank(
  den_collection::FeynmanDenominatorCollection
)::Int
  loop_momenta = map(Basic, den_collection.loop_momenta)
  independent_external_momenta = map(Basic, den_collection.external_momenta[begin:end-1])
  den_list = map(Basic, den_collection.denominators)
  return calculate_denominator_collection_rank(
    loop_momenta, independent_external_momenta, den_list
  )
end # function calculate_denominator_collection_rank
function calculate_denominator_collection_rank(
  loop_momenta::Vector{Basic},
  independent_external_momenta::Vector{Basic},
  den_list::Vector{FeynmanDenominator},
)::Int
  n_loop = length(loop_momenta)
  n_ind_ext = length(independent_external_momenta)
  n_sp = binomial(n_loop, 2) + n_loop * (n_ind_ext + 1)
  # @show n_sp

  denominator_exprs = map((expand ∘ Basic), get_denominator_expr(den_list))

  # all_independent_scalar_products = Vector{Basic}(undef, n_sp)

  coeff_mat = Matrix{Basic}(undef, length(denominator_exprs), n_sp)
  for (ii, expr) ∈ enumerate(denominator_exprs)
    jj = 1
    for ii_loop ∈ 1:n_loop
      tmp_coeff = SymEngine.coeff(expr, loop_momenta[ii_loop], Basic(1))

      for jj_loop ∈ ii_loop:n_loop
        coeff_mat[ii, jj] = if ii_loop == jj_loop
          SymEngine.coeff(expr, loop_momenta[ii_loop], Basic(2))
        else
          SymEngine.coeff(tmp_coeff, loop_momenta[jj_loop], Basic(1))
        end
        jj += 1
      end # for jj_loop

      for jj_ext ∈ 1:n_ind_ext
        coeff_mat[ii, jj] = SymEngine.coeff(tmp_coeff, independent_external_momenta[jj_ext], Basic(1))
        jj += 1
      end # for jj_ext
    end # for ii_loop
    @assert jj == n_sp + 1
  end # for (ii, expr)
  # @show coeff_mat
  return (rank ∘ map)(N, coeff_mat)
end # function calculate_denominator_collection_rank

###########################################
# greedy algorithm
function get_greedy_cover(
  original_indices_list::Vector{Vector{Int}}
)::Vector{Vector{Int}}
###########################################
  # length(indices_list) == 1 && return indices_list
  universe = reduce(union, original_indices_list; init=Int[])
  universe_copy = copy(universe)
  indices_list = deepcopy(original_indices_list)

  new_indices_list = Vector{Int}[]

  while !isempty(universe)
    _, pos = findmax(indices -> length(indices ∩ universe), indices_list)
    # @show pos
    push!(new_indices_list, copy(indices_list[pos]))
    # @show new_indices_list
    setdiff!(universe, indices_list[pos])
    deleteat!(indices_list, pos)
  end # while

  @assert (isempty ∘ symdiff)(universe_copy, reduce(union, new_indices_list))

  # @show new_indices_list
  return new_indices_list
end # function greedy

function write_topology(
  topology_dict::Dict{FeynmanDenominatorCollection, Dict{Int, Dict{Basic, Basic}}},
  amp_file_list::Vector{String}, topology_directory::String, kinematic_relations::Dict{Basic, Basic}
)::Nothing
  counter = 1
  for (key_collection, covering_dict) ∈ topology_dict
    topology_filename = joinpath(topology_directory, "topology$(counter).jld2")
    @info "Writing topology file: $topology_filename."
    jldopen(topology_filename, "w+") do f
      f["loop_momenta"] = key_collection.loop_momenta
      f["external_momenta"] = key_collection.external_momenta
      f["kinematic_relation"] = Dict{String, String}(
        string(key) => string(value) for (key, value) ∈ kinematic_relations
      )
      f["covering_amplitudes"] = Dict{String, Dict{String, String}}(
        relative_path(topology_filename, amp_file_list[key]) => Dict{String, String}(string(k) => string(v) for (k, v) ∈ value)
          for (key, value) ∈ covering_dict
      )
      f["denominators"] = map(string, key_collection.denominators)
    end # jldopen

    topology_outname = joinpath(topology_directory, "topology$(counter).out")
    open(topology_outname, "w") do f
      write(f, key_collection)
      write(f, "    covering amplitude files:\n")
      for (key, value) ∈ covering_dict
        write(f, "        $(relative_path(topology_outname, amp_file_list[key])):")
        isempty(value) && write(f, " No momentum shifts required.")
        write(f, "\n")
        for (k, v) ∈ value
          write(f, "            - $k => $v\n")
        end # for
      end # for
    end # open

    counter += 1
  end # for
end # function write_topology
