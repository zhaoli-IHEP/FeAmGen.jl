function construct_topology(
  amp_dir::String;
  method::Symbol=:Canonicalization,
  options::Dict{String, <:Any}=Dict{String, Any}()
)
  # call old construct_den_topology API #######################################
  if method == :Canonicalization
    mom_shift_opt = haskey(options, "mom_shift_opt") ?
      options["mom_shift_opt"] : false
    mom_shift_opt = isa(mom_shift_opt, Bool) ?
      mom_shift_opt : false
    ref_dentop_collect = haskey(options, "ref_dentop_collect") ?
      options["ref_dentop_collect"] : DenTop[]
    ref_dentop_collect = isa(
        ref_dentop_collect, Union{String,Vector{DenTop}}
      ) ? ref_dentop_collect : DenTop[]
    return construct_den_topology(amp_dir;
      mom_shift_opt=mom_shift_opt,
      ref_dentop_collect=ref_dentop_collect
    )
  end # if
  #############################################################################

  return construct_topology(Val(method), amp_dir, options)
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
  original_den_collection_list = read_loop_denominators(Val(:Amplitude), amp_filename_list)

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

read_loop_denominators(opt::Any, file_name_list::Vector{String}) =
  map(file_name -> read_loop_denominators(opt, file_name), file_name_list)
function read_loop_denominators(::Val{:Amplitude}, amplitude_file_name::String)::FeynmanDenominatorCollection
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

  return FeynmanDenominatorCollection(loop_momenta, external_momenta, den_list)
end # function read_loop_denominators

function write_topology(
  topology_dict::Dict{FeynmanDenominatorCollection, Dict{Int, Vector{Dict{Basic, Basic}}}},
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
      f["covering_amplitudes"] = Dict{String, Vector{Dict{String, String}}}(
        relative_path(topology_filename, amp_file_list[key]) => Dict{String, String}[
          Dict{String, String}(string(k) => string(v) for (k, v) ∈ repl_rules)
            for repl_rules ∈ value
        ] for (key, value) ∈ covering_dict
      )
      f["denominators"] = map(string, key_collection.denominators)
    end # jldopen

    topology_outname = joinpath(topology_directory, "topology$(counter).out")
    open(topology_outname, "w") do f
      write(f, key_collection)
      write(f, "    covering amplitude files:\n")
      for (key, value) ∈ covering_dict
        write(f, "        $(relative_path(topology_outname, amp_file_list[key])):\n")
        for repl_rules ∈ value
          write(f, "            - |\n")
          for (k, v) ∈ repl_rules
            write(f, "                - $k => $v\n")
          end # for
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
