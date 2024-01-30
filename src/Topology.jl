include("Topology_legacy.jl")
include("Topology_Canon.jl")
include("Topology_Pak.jl")

function construct_den_topology(
  amp_dir::String;
  method_to_find_momentum_shifts::Symbol=:Canonicalization, # :Canonicalization or :PakAlgorithm
  check_momentum_shifts::Bool=true,
  use_reference_topology::Bool=false,
  reference_topology_directory::String="",
  reference_topologies::Vector{FeynmanDenominatorCollection}=FeynmanDenominatorCollection[],
  check_all_kinematic_relations::Bool=false,
  find_external_momentum_shifts::Bool=false,
  recheck_momentum_shifts::Bool=false
)::Vector{FeynmanDenominatorCollection}

  # check options #############################################################
  @assert method_to_find_momentum_shifts ∈ [:Canonicalization, :PakAlgorithm] "Do not support method $method_to_find_momentum_shifts to find momentum shifts."

  use_reference_topology &&
    return construct_den_topology_with_reference_topology(
      amp_dir;
      method_to_find_momentum_shifts=method_to_find_momentum_shifts,
      check_momentum_shifts=check_momentum_shifts,
      reference_topology_directory=reference_topology_directory,
      reference_topologies=reference_topologies,
      check_all_kin_relation=check_all_kin_relation,
      find_external_momentum_shifts=find_external_momentum_shifts
    )
  # end check options #########################################################

  # init ######################################################################
  amp_filename_list, ds_list, kin_relation_dict  = read_loop_denominators(
    Val(:AmplitudeDirectory), amp_dir, check_all_kinematic_relations
  )
  # end init ##################################################################

  # step 1 ####################################################################
  # construct a minimal set of Feynman integral topologies ####################
  @info "Minimalizing the initial topologies."
  minimal_topology_list, covering_indices = minimize_topology_list_directly(
    Val(method_to_find_momentum_shifts), ds_list
  ) # end minimal_topology_list
  @info "Done and got $(length(minimal_topology_list)) initial topologies."
  # end step 1 ################################################################

  # step 2 ####################################################################
  result_topology_list, result_repl_rules_list = if method_to_find_momentum_shifts == :Canonicalization
    nothing # TBD
  elseif method_to_find_momentum_shifts == :PakAlgorithm
    construct_den_topology(
      Val(:PakAlgorithm), minimal_topology_list,
      kin_relation_dict,
      find_external_momentum_shifts
    ) # end construct_den_topology
  end # if
  covering_dict = Dict{FeynmanDenominatorCollection, Dict{Int, Dict{Basic, Basic}}}()
  for ii ∈ eachindex(result_topology_list)
    amp_repl_rules_dict = Dict{Int, Dict{Basic, Basic}}()
    result_repl_rules = result_repl_rules_list[ii]
    for jj ∈ keys(result_repl_rules), kk ∈ covering_indices[jj]
      amp_repl_rules_dict[kk] = result_repl_rules[jj]
    end # for kk
    covering_dict[result_topology_list[ii]] = amp_repl_rules_dict
  end # for ii
  # end step 2 ################################################################

  # step 3 ####################################################################
  # complete topologies #######################################################
  if method_to_find_momentum_shifts == :Canonicalization
    nothing # TBD
  elseif method_to_find_momentum_shifts == :PakAlgorithm
    nothing # TBD
  end # if
  # end step 3 ################################################################

  # recheck momentum shifts #####################################################
  if recheck_momentum_shifts
    @info "Re-checking all original topologies are covered by the constructed topologies."
    covering_dict = Dict{FeynmanDenominatorCollection, Dict{Int, Dict{Basic, Basic}}}()
    for (ii, topology) ∈ enumerate(result_topology_list)
      repl_rules = Dict{Int, Dict{Basic, Basic}}()
      for (jj, den_collection) ∈ enumerate(ds_list)
        @info "Checking if the constructed topology #$ii (total: $(length(result_topology_list))) covering the original topology #$jj (total: $(length(ds_list)))."
        if method_to_find_momentum_shifts == :Canonicalization
          nothing # TBD
        elseif method_to_find_momentum_shifts == :PakAlgorithm
          momentum_shifts = find_Pak_momentum_shifts(
            topology, den_collection;
            SP_replacements=kin_relation_dict,
            find_external_momentum_shifts=find_external_momentum_shifts
          )
          isnothing(momentum_shifts) && continue
          repl_rules[jj] = momentum_shifts
        end # if
      end # for (jj, den_collection)
      covering_dict[topology] = repl_rules
    end # for topology

    @assert (isempty ∘ symdiff)(
      reduce(union, map(collect ∘ keys, (collect ∘ values)(covering_dict))),
      eachindex(ds_list)
    )
  end # if
  # end recheck momentum shifts ###############################################


  # write topology ############################################################
  top_dir_list = splitpath(amp_dir)
  topology_name = if method_to_find_momentum_shifts == :Canonicalization
    "Canonicalization_topologies"
  elseif method_to_find_momentum_shifts == :PakAlgorithm
    "Pak_topologies"
  end # if
  top_dir_list[end] = if endswith(top_dir_list[end], "amplitudes")
    replace(top_dir_list[end], "amplitudes" => topology_name)
  else
    top_dir_list[end] *= "_" * topology_name
  end
  topology_directory = joinpath(top_dir_list)
  bk_mkdir(topology_directory)
  write_topology(covering_dict, amp_filename_list, topology_directory, kin_relation_dict)
  # end write topology ########################################################

  return result_topology_list
end # function construct_topology

function construct_topology(
  ::Val{:PakAlgorithm}, amp_dir::String,
  options::Dict{String, <:Any}=Dict{String, Any}()
)::Dict
  # load den_collection #######################################################
  @assert isdir(amp_dir) "Got non-existent directory: $amp_dir."
  amp_filename_list = filter!(endswith(".jld2"), readdir(amp_dir; join=true, sort=false))
  sort!(amp_filename_list, by=filename -> begin
                                            main_name, _ = (splitext ∘ basename)(filename)
                                            index_str = match(r"[1-9]\d*", main_name).match
                                            parse(Int, index_str)
                                          end)
  @assert !isempty(amp_filename_list) "Cannot find any amplitude file with extension .jld2 in $amp_dir."
  original_den_collection_list = read_loop_denominators(Val(:AmplitudeFile), amp_filename_list)

  # load kinematic relation ###################################################
  kin_relation_dict = Dict{Basic, Basic}()
  if haskey(options, "check_all_kin_relation") && options["check_all_kin_relation"]
    for amp_file ∈ amp_filename_list
      kin_relation = load(amp_file, "kin_relation")
      for (key, value) ∈ kin_relation
        key = Basic(key)
        value = Basic(value)
        get_name(key) == "SP" || continue
        haskey(kin_relation_dict, key) && 
          @assert kin_relation_dict[key] == value """
          \rGot different kinematic relation:
          \r    - $key => $(kin_relation_dict[key]) in previous files;
          \r    - $key => $value in $amp_file.
          """
        kin_relation_dict[key] = value
      end # for (key, value)
    end # for amp_file
  end # if
  for (key, value) ∈ load(first(amp_filename_list), "kin_relation")
    key = Basic(key)
    get_name(key) == "SP" || continue
    kin_relation_dict[key] = Basic(value)
  end # for (key, value)

  # construct topology ########################################################
  covering_dict = Dict{FeynmanDenominatorCollection, Dict{Int, Vector{Dict{Basic, Basic}}}}()
  remaining_indices = (collect ∘ eachindex)(original_den_collection_list)
  while !isempty(remaining_indices)
    @info "Number of remaining indices: $(length(remaining_indices))."

    _, index_index = findmax(length, original_den_collection_list[remaining_indices])
    largest_index = remaining_indices[index_index]
    # index_index, largest_index = 1, first(remaining_indices)
    while true
      tmp_den_collection = original_den_collection_list[largest_index]
      tmp_index_index = findnext(
        index -> tmp_den_collection ≠ original_den_collection_list[index] &&
          !is_Pak_equivalent(
            tmp_den_collection, original_den_collection_list[index];
            SP_replacements=kin_relation_dict
          ) && original_den_collection_list[index] ⊇ tmp_den_collection || check_Pak_covering(
            original_den_collection_list[index], tmp_den_collection;
            SP_replacements=kin_relation_dict
          ),
        remaining_indices, index_index + 1
      )
      isnothing(tmp_index_index) && break
      index_index, largest_index = tmp_index_index, remaining_indices[tmp_index_index]
    end

    key_collection = original_den_collection_list[largest_index]
    all_repl_rules_list = Vector{Vector{Dict{Basic, Basic}}}(undef, length(remaining_indices))
    for (ii, index) ∈ enumerate(remaining_indices)
      target_collection = original_den_collection_list[index]
      if index == largest_index || key_collection ⊇ target_collection
        all_repl_rules_list[ii] = [Dict{Basic, Basic}()] # do not check the same topology or directly sub-topology
        continue
      end # if
      print("\rChecking $largest_index covering: $(ii) / $(length(remaining_indices)).")
      all_repl_rules_list[ii] = find_Pak_momentum_shifts(
        key_collection, target_collection;
        SP_replacements=kin_relation_dict,
        find_external_momentum_shifts_flag=haskey(options, "find_external_momentum_shifts") ?
          options["find_external_momentum_shifts"] : false
      )
    end # for index ∈ remaining_indices
      
    indices = findall(!isempty, all_repl_rules_list)
    covering_dict[key_collection] = Dict{Int, Vector{Dict{Basic, Basic}}}(
      remaining_indices[indices] .=> all_repl_rules_list[indices]
    )
    deleteat!(remaining_indices, indices)
    print("\r")
  end # while

  top_dir_list = splitpath(amp_dir)
  top_dir_list[end] = if endswith(top_dir_list[end], "amplitudes")
    replace(top_dir_list[end], "amplitudes" => "Pak_topologies")
  else
    top_dir_list[end] *= "_Pak_topologies"
  end
  topology_directory = joinpath(top_dir_list)
  bk_mkdir(topology_directory)
  write_topology(covering_dict, amp_filename_list, topology_directory, kin_relation_dict)

  return covering_dict
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

  original_den_list_list = map(dc -> unique(dc.denominators), den_collection_list)
  # end initial information ###################################################

  den_list_list = Vector{FeynmanDenominator}[]
  while !isempty(original_den_list_list)
    _, index = findmax(length, original_den_list_list)
    target_den_list = original_den_list_list[index]
    push!(den_list_list, target_den_list)
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
        isempty(value) && write(f, "No momentum shifts required.")
        write(f, "\n")
        for (k, v) ∈ value
          write(f, "            - $k => $v\n")
        end # for
      end # for
    end # open

    counter += 1
  end # for
end # function write_topology

function relative_path(base_path::String, target_path::String)::String
  @assert ispath(base_path) "Got non-existent path: $base_path."
  @assert ispath(target_path) "Got non-existent path: $target_path."

  base_abspath = abspath(base_path)
  target_abspath = abspath(target_path)
  if base_abspath == target_abspath
    isfile(base_abspath) && return joinpath(".", basename(base_abspath))
    isdir(base_abspath) && return "."
  end # if

  base_abspath_list = splitpath(base_abspath)
  target_abspath_list = splitpath(target_abspath)

  while true
    isempty(base_abspath_list) && return joinpath(".", target_abspath_list...)
    isempty(target_abspath_list) && break
    first(base_abspath_list) == first(target_abspath_list) || break
    popfirst!(base_abspath_list)
    popfirst!(target_abspath_list)
  end # while

  result_path_list = if isfile(base_abspath)
    vcat([".." for _ ∈ base_abspath_list[begin:end-1]], target_abspath_list)
  else
    @assert isdir(base_abspath)
    vcat([".." for _ ∈ base_abspath_list], target_abspath_list)
  end # if

  return joinpath(result_path_list)
end # function relative_path
