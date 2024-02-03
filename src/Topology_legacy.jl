import Base: isempty, issubset, length, union

###########################################
struct DenTop
  n_loop::Int
  indep_ext_mom::Vector{Basic}
  den_list::Vector{Basic}
end # mutable struct DenTop
###########################################




isempty( dentop::DenTop ) = isempty( dentop.den_list )

issubset( dentop1::DenTop, dentop2::DenTop )::Bool =
    dentop1.n_loop ≤ dentop2.n_loop &&
    issubset( dentop1.indep_ext_mom, dentop2.indep_ext_mom ) &&
    issubset( dentop1.den_list, dentop2.den_list )

length( dentop::DenTop )::Int64 = length( dentop.den_list )

###########################################
function union(
    dentop1::DenTop,
    dentop2::DenTop
)::DenTop
###########################################
  @assert dentop1.n_loop == dentop2.n_loop
  @assert dentop1.indep_ext_mom == dentop2.indep_ext_mom

  new_den_list = union( dentop1.den_list, dentop2.den_list )

  return DenTop( dentop1.n_loop, dentop1.indep_ext_mom, new_den_list )
end # function union

###########################################
function union(
    dentop::DenTop,
    den_list::Vector{Basic}
)::DenTop
###########################################

  @assert (get_loop_index∘last∘get_loop_momenta)( den_list ) ≤ dentop.n_loop
  @assert get_ext_momenta( den_list ) ⊆ dentop.indep_ext_mom

  return DenTop( dentop.n_loop, dentop.indep_ext_mom, union(dentop.den_list,den_list) )
end # function union

union( dentop::DenTop, den::Basic ) = union( dentop, [den] )

union( dentop::DenTop )::DenTop =
  DenTop( dentop.n_loop, dentop.indep_ext_mom, union(dentop.den_list) )

union( dentop::DenTop, dentop_list... )::DenTop =
  union( dentop, reduce( union, dentop_list ) )


###########################################
function is_valid_dentop(
    dentop::DenTop
)::Bool
###########################################

  n_loop = dentop.n_loop
  indep_ext_mom = dentop.indep_ext_mom
  den_list = dentop.den_list

  !all( is_ext_mom, indep_ext_mom ) && return false
  if !isempty(den_list)
    unique_den_list = unique( den_list )
    !all( is_FunctionSymbol, unique_den_list ) && return false
    !all( den->get_name(den)=="Den", unique_den_list ) && return false

    loop_momenta = get_loop_momenta( unique_den_list )
    ext_momenta = get_ext_momenta( unique_den_list )

    max_loop_index = isempty( loop_momenta ) ? 0 : (get_loop_index∘last∘get_loop_momenta)( unique_den_list )
    max_ext_index = isempty( ext_momenta ) ? 0 : (get_ext_index∘last∘get_ext_momenta)( unique_den_list )
    max_ind_ext_index = isempty( indep_ext_mom ) ? 0 : (first∘findmax∘map)( get_ext_index, indep_ext_mom )

    max_loop_index > n_loop && return false
    max_ext_index > max_ind_ext_index && return false
  end # if

  return true
end # function is_valid_dentop














###########################################
function gen_sp_dict(
  dentop::DenTop
)::Dict{Basic, Basic}
###########################################

  n_loop = dentop.n_loop
  n_ext_mom = length(dentop.indep_ext_mom)
  sp_index = 1
  sp_dict = Dict{Basic, Basic}()

  for ii ∈ 1:n_loop, jj ∈ ii:dentop.n_loop
    qi, qj = Basic("q$ii"), Basic("q$jj")
    sp_dict[ make_SP(qi, qj) ] = Basic("sp$(sp_index)")
    sp_index += 1
  end # for ii, jj

  for one_ext_mom ∈ dentop.indep_ext_mom, loop_ii ∈ 1:n_loop
    q = Basic("q$(loop_ii)")
    sp_dict[ make_SP(one_ext_mom, q) ] = Basic("sp$(sp_index)")
    sp_index += 1
  end # one_ext_mom, q

  @assert sp_index == n_loop * (n_loop + 1) / 2 + n_loop * n_ext_mom + 1

  return sp_dict

end # function gen_sp_dict

###########################################
function get_vac_den_list(
    dentop::DenTop
)::Vector{Basic}
###########################################

  den_list = dentop.den_list
  ext_momenta = get_ext_momenta( den_list )
  vac_den_list = subs.( den_list, (ext_momenta .=> 0)... )
  unique!(abs, vac_den_list)

  return vac_den_list

end # function get_vac_den_list

###########################################
function get_coeff_mat_mom2_sp(
    dentop::DenTop
)::Matrix{Rational}
###########################################

  sp_dict = gen_sp_dict( dentop )
  n_sp = length(sp_dict)

  mom2_list = subs.( map( make_SP∘expand∘first∘get_args, dentop.den_list ), Ref(sp_dict) )

  return map( Rational∘Int,
              coefficient_matrix( mom2_list, [ Basic("sp$(index)") for index ∈ 1:n_sp ] ) )

end # function get_coeff_mat_mom2_sp

###########################################
function get_superior_dentop_collect(
  dentop_collect::Vector{DenTop}
)::Vector{DenTop}
###########################################

  new_dentop_collect = DenTop[]

  for dentop ∈ dentop_collect
    included_by_pos = findfirst( new_dentop->dentop⊆new_dentop, new_dentop_collect )
    !isnothing(included_by_pos) && continue
    filter!( new_dentop->new_dentop⊈dentop, new_dentop_collect )
    push!( new_dentop_collect, dentop )
  end # for dentop

  return new_dentop_collect

end # function get_superior_dentop_collect

###########################################
function get_cover_indices_list(
    dentop_collect::Vector{DenTop}
)::Vector{Vector{Int}}
###########################################
  n_loop = first( dentop_collect ).n_loop
  n_ind_ext = length( first(dentop_collect).indep_ext_mom )
  n_sp::Int = (n_loop + 1) * n_loop / 2 + n_loop * n_ind_ext
  preferred_vac_mom_list = if haskey( preferred_vac_mom_dict(), n_loop )
    preferred_vac_mom_dict()[n_loop]
  else
    Vector{Basic}[]
  end

  n_dentop = length( dentop_collect )
  graded_indices_list = [ [ [ii] for ii ∈ 1:n_dentop ] ]
  prev_indices_list = last( graded_indices_list )

  for _ ∈ 2:n_dentop
    this_indices_list = Vector{Int}[]
    for one_indices ∈ prev_indices_list
      for one_index ∈ (last(one_indices)+1):n_dentop
        push!( this_indices_list, (sort∘union)( one_indices, one_index ) )
      end # for one_index
    end # for one_indices

    dentop_union_list = map( indices->union(dentop_collect[indices]...),
                              this_indices_list )

    pos_list = findall( dentop->begin
                                  vac_den_list = get_vac_den_list(dentop)
                                  vac_mom_list = map( first∘get_args, vac_den_list )
                                  flag = length(dentop) ≤ n_sp && length(vac_den_list) ≤ 3 * (n_loop - 1)
                                  flag &= isempty( preferred_vac_mom_list ) ? true :
                                            any( mom_list->vac_mom_list⊆mom_list, preferred_vac_mom_list )
                                  # flag &= length(dentop) == (rank∘get_coeff_mat_mom2_sp)(dentop)
                                end, # begin
                          dentop_union_list )
    this_indices_list = this_indices_list[ pos_list ]
    dentop_union_list = dentop_union_list[ pos_list ]

    pos_list = findall( dentop->length(dentop)==(rank∘get_coeff_mat_mom2_sp)(dentop),
                          dentop_union_list )

    isempty( pos_list ) && break

    prev_indices_list = this_indices_list[ pos_list ]
    push!( graded_indices_list, prev_indices_list )
  end # for

  result_indices_list = Vector{Int}[]
  for n_index ∈ (reverse∘eachindex)( graded_indices_list )
    this_indices_list = graded_indices_list[ n_index ]
    filter!( this_indices->all(result_indices->isempty(result_indices∩this_indices),
                            result_indices_list),
                this_indices_list )
    isempty( this_indices_list ) && continue
    @assert all( indices->length(indices)==n_index, this_indices_list )
    append!( result_indices_list, greedy( this_indices_list ) )
  end

  return result_indices_list

end # function get_cover_indices_list

###########################################
# greedy algorithm
function greedy(
    indices_list::Vector{Vector{Int}}
)::Vector{Vector{Int}}
###########################################
  universe = reduce(union,  indices_list)

  new_indices_list = Vector{Int}[]

  while !isempty(universe)
    _, pos = findmax( indices->length(indices∩universe), indices_list )
    push!( new_indices_list, indices_list[pos] )
    setdiff!( universe, indices_list[pos] )
    deleteat!( indices_list, pos )
  end # while

  return new_indices_list
end # function greedy

###########################################
function make_complete_dentop_collect(
  dentop_list::Vector{DenTop}
)::Vector{DenTop}
###########################################

  n_loop = first(dentop_list).n_loop
  @assert haskey( preferred_vac_mom_dict(), n_loop ) "$(n_loop)-loop is not supported now."

  indep_ext_mom = first(dentop_list).indep_ext_mom
  n_sp::Int = (n_loop + 1) * n_loop / 2 + n_loop * length( indep_ext_mom )
  preferred_vac_mom_lists = preferred_vac_mom_dict()[n_loop]

  incomplete_dentop_list = copy(dentop_list)
  complete_dentop_list = DenTop[]

  cover_indices_list = get_cover_indices_list( incomplete_dentop_list )
  @assert (length∘reduce)( union, cover_indices_list ) == length( incomplete_dentop_list )

  for indices ∈ cover_indices_list
    to_be_complete_dentop = union( incomplete_dentop_list[indices]... )

    while length(to_be_complete_dentop) < n_sp
      # println( "$to_be_complete_dentop need to be completed." )
      # @assert length(to_be_complete_dentop) == (rank∘get_coeff_mat_mom2_sp)( to_be_complete_dentop )
      vac_mom_list = begin
        this_vac_den_list = get_vac_den_list(to_be_complete_dentop)
        this_vac_mom_list = map( first∘get_args, this_vac_den_list )
        selected_index = findfirst( vac_mom_list->this_vac_mom_list⊆vac_mom_list,
                                      preferred_vac_mom_lists )
        @assert !isnothing(selected_index)
        preferred_vac_mom_lists[ selected_index ]
      end # vac_mom_list

      for ext_mom ∈ vcat( zero(Basic), indep_ext_mom ), q ∈ vac_mom_list, the_sign ∈ [1,-1]
        trial_den = SymFunction("Den")(expand(q + the_sign * ext_mom), 0, 0)
        trial_top = union( to_be_complete_dentop, trial_den )
        rank_trial_top = (rank∘get_coeff_mat_mom2_sp)(trial_top)
        if rank_trial_top > length(to_be_complete_dentop)
          to_be_complete_dentop = trial_top
          rank_trial_top == n_sp && break
        end # if
      end # for ext_mom
    end # while

    @assert length( to_be_complete_dentop ) == (rank∘get_coeff_mat_mom2_sp)( to_be_complete_dentop ) == n_sp
    # println( "$to_be_complete_dentop is complete." )
    # println()
    push!( complete_dentop_list, to_be_complete_dentop )
  end # for indices

  complete_dentop_list = get_superior_dentop_collect( complete_dentop_list )

  return complete_dentop_list

end # function make_complete_dentop_collect

###########################################
function readin_dentop_collect( dentop_dir::String )::Vector{DenTop}
###########################################
  @assert isdir( dentop_dir )
  dentop_file_list = readdir( dentop_dir; join=true, sort=false )
  filter!( isfile, dentop_file_list )
  filter!( contains(r"^topology[1-9]\d*.jld2$") ∘ basename, dentop_file_list )
  sort!( dentop_file_list; by=get_diagram_index )

  dentop_collect = DenTop[]
  for dentop_file ∈ dentop_file_list
    dentop = load( dentop_file )
    push!( dentop_collect, DenTop( dentop["n_loop"],
                                    to_Basic( dentop["indep_ext_mom"] ),
                                    to_Basic( dentop["den_list"] ) ) )
  end # for dentop_file

  return dentop_collect

end # function readin_dentop_collect

###########################################
function find_mom_shift(
    den_list::Vector{Basic},
    ref_dentop_collect::Vector{DenTop}
)::Dict{Basic,Basic}
###########################################
  isempty( ref_dentop_collect ) && return Dict{Basic,Basic}()
  @assert all( den->get_name(den)=="Den"&&(length∘get_args)(den)==3, den_list)
  @funs Den

  # get information
  tmp_q_list = get_loop_momenta( den_list )
  k_list = get_ext_momenta( den_list )
  n_loop = isempty(tmp_q_list) ? 0 : (get_loop_index∘last)( tmp_q_list )
  preffered_vac_mom_lists = preferred_vac_mom_dict()[n_loop]
  q_list = [ Basic("q$ii") for ii ∈ 1:n_loop ]
  diff_q_list = setdiff(q_list, tmp_q_list)
  !isempty(diff_q_list) && @warn """
  Add $(join(diff_q_list, ", ")) to the loop momenta list.
  """
  q_null_dict = Dict( q_list .=> zero(Basic) )

  # find all branches of the loop momenta.
  mom_list = map( normalize_loop_mom_single∘first∘get_args,
                  den_list )
  mom_q_coeff_mat = coefficient_matrix( mom_list, q_list )
  @assert (isempty∘setdiff)( unique(mom_q_coeff_mat), Basic[-1,0,1] )
  mom_q_coeff_list = (unique∘collect∘eachrow)( mom_q_coeff_mat )
  sort!( mom_q_coeff_list; 
          by=coeff_list->begin
            num_q = (length∘findall)(!iszero,coeff_list)
            new_coeff_list = map( coeff->subs(coeff,Dict(Basic(-1)=>2,Basic(0)=>3)), coeff_list )
            num_q, new_coeff_list
          end )
  vac_mom_list = [ sum( mom_q_coeff .* q_list ) for mom_q_coeff ∈ mom_q_coeff_list ]
  branch_mom_indices_list = Vector{Int}[]
  for mom_q_coeff ∈ mom_q_coeff_list
    mom_indices = findall( ==(mom_q_coeff), eachrow(mom_q_coeff_mat) )
    push!( branch_mom_indices_list, mom_indices )
  end # for mom_q_coeff

  for selected_branch ∈ permutations( eachindex(vac_mom_list), length(q_list) )
    for sign_list ∈ Iterators.product( [ (1, -1) for _ ∈ q_list ]... )

      # check invertible q coeff matrix
      selected_vac_mom_list = vac_mom_list[ selected_branch ]
      selected_coeff_mat = coefficient_matrix( selected_vac_mom_list, q_list )
      (iszero∘expand∘get_det)( selected_coeff_mat ) && break
      inv_selected_coeff_mat = inv( selected_coeff_mat )

      # check vacumm momenta is preferred
      vac_repl_rule = Dict( q_list .=> inv_selected_coeff_mat * ( q_list .* sign_list ) )
      new_vac_mom_list = map( mom->(expand∘subs)(mom,vac_repl_rule), vac_mom_list )
      map!( normalize_loop_mom_single, new_vac_mom_list, new_vac_mom_list )
      !any( preffered_vac_mom_list->new_vac_mom_list⊆preffered_vac_mom_list,
            preffered_vac_mom_lists ) && continue
      # @info "checked vacumm momenta is preferred."

      # construct the repl rule
      for selected_mom_indices ∈ Iterators.product( branch_mom_indices_list[selected_branch]... )
        selected_ext_mom_list = map( mom->subs(mom,q_null_dict), mom_list[collect(selected_mom_indices)] )
        repl_rule = Dict( q_list .=> map( expand, inv_selected_coeff_mat * (sign_list .* q_list .- selected_ext_mom_list) ) )
        new_mom_list = map( mom->(expand∘subs)(mom,repl_rule), mom_list )
        map!( normalize_loop_mom_single, new_mom_list, new_mom_list )
        (!isempty∘setdiff)( (unique∘map)( abs, coefficient_matrix( new_mom_list, k_list ) ),
                              Basic[0,1] ) && continue

        # check cover
        new_den_list = [  begin
                            _, mass, ieta = get_args(den)
                            Den(new_mom, mass, ieta)
                          end for (new_mom, den) ∈ zip( new_mom_list, den_list ) ]
        new_dentop = DenTop( n_loop, k_list, new_den_list )
        if any( dentop->new_dentop⊆dentop||dentop⊆new_dentop, ref_dentop_collect )
          for key ∈ (collect∘keys)( repl_rule )
            repl_rule[key] == key && delete!( repl_rule, key )
          end # for key
          return repl_rule
        end # if
      end # for selected_mom_indices
    end # for sign_list
  end # for selected_branch_indices_list

  return Dict{Basic,Basic}()
end # function find_mom_shift

###########################################
function construct_den_topology_legacy(
    amp_dir::String; 
    mom_shift_opt::Bool=false,
    ref_dentop_collect::Union{String,Vector{DenTop}}=DenTop[],
)::Vector{DenTop}
###########################################

  @assert isdir(amp_dir)
  endswith( amp_dir, '/') && (amp_dir = amp_dir[begin:end-1])
  @assert endswith( amp_dir, "_amplitudes" )
  root_dir, shifted_amp_dir, topology_dir = begin
    root_dir, amp_name = (splitdir∘abspath)(amp_dir)

    root_dir,
      joinpath( root_dir, replace( amp_name, "_amplitudes" => "_shifted_amplitudes" ) ),
      joinpath( root_dir, replace( amp_name, "_amplitudes" => mom_shift_opt ? "_shifted_topologies" : "_topologies" ) )
  end # topology_dir
  mom_shift_opt && bk_mkdir( shifted_amp_dir )
  bk_mkdir( topology_dir )
  dentop_collect = !mom_shift_opt ? DenTop[] :
                    ( isa(ref_dentop_collect, String) ?
                        readin_dentop_collect(ref_dentop_collect) :
                        ref_dentop_collect )
  @info "Read in $(length(dentop_collect)) reference topologies."

  amp_file_list = readdir( amp_dir; join=true, sort=false )
  filter!( isfile, amp_file_list )
  filter!( contains(r"^amp[1-9]\d*.jld2$") ∘ basename, amp_file_list )
  amp_file_indices = map( get_diagram_index, amp_file_list )
  amp_order = sortperm( amp_file_indices )
  amp_file_list = amp_file_list[amp_order]
  shifted_amp_file_list = mom_shift_opt ? map( file_name->joinpath(shifted_amp_dir,replace(basename(file_name),r"^amp"=>"shifted_amp")),
                                                amp_file_list ) : String[]
  amp_file_indices = amp_file_indices[amp_order]
  @info "Found $(length(amp_file_list)) amplitude files at $amp_dir."

  n_loop, ext_mom_list = jldopen( first(amp_file_list), "r" ) do jld_file
    jld_file["n_loop"], to_Basic( jld_file["ext_mom_list"] )
  end # n_loop, ext_mom_list
  @assert all( amp_file->jldopen( amp_file, "r" ) do jld_file
                            jld_file["n_loop"] == n_loop &&
                            to_Basic( jld_file["ext_mom_list"] ) == ext_mom_list
                          end, amp_file_list[2:end] ) "The amplitude files are not consistent!"
  if iszero(n_loop)
    @warn "There is nothing to do for tree-level!"
    return DenTop[]
  end # if
  indep_ext_mom = ext_mom_list[begin:end-1]
  mom_shift_opt && @assert all( ref_dentop->ref_dentop.n_loop==n_loop&&ref_dentop.indep_ext_mom==indep_ext_mom,
                                dentop_collect ) "The reference den-topologies are not consistent with the amplitudes!"

  #--------------------------------
  jld_file = jldopen( first(amp_file_list), "r" ) 
  mom_symmetry = to_Basic_dict( jld_file["mom_symmetry"] ) 
  indep_ext_mom = (unique∘free_symbols∘map)( x->subs(x,mom_symmetry), indep_ext_mom ) 
  close( jld_file )
  #--------------------------------

  backup_dentop_collect = DenTop[]
  mom_shift_collect = Dict{Basic,Basic}[]
  for (index, amp_index) ∈ enumerate( amp_file_indices )
    amp_file = amp_file_list[index]
    shifted_amp_file = mom_shift_opt ? shifted_amp_file_list[index] : 0

    jld_file = jldopen( amp_file, "r" )

    den_list = to_Basic( jld_file["loop_den_list"] )

    if mom_shift_opt
      @info "Finding momenta shift for $amp_index/$(last(amp_file_indices))..."

      amp_list = to_Basic( jld_file["amp_lorentz_list"] )
      # mom_list = map( first∘get_args, den_list )

      mom_shift = find_mom_shift( den_list, dentop_collect )
      map!( den->subs(den,mom_shift), den_list, den_list )
      map!( amp->subs(amp,mom_shift), amp_list, amp_list )
      den_list = normalize_loop_mom( den_list )
      # mom_list = map( first∘get_args, den_list )
      dentop = DenTop( n_loop, indep_ext_mom, den_list )
      push!( backup_dentop_collect, dentop )
      push!( mom_shift_collect, mom_shift )

      covered_indices = findfirst( ref_dentop->dentop⊆ref_dentop, dentop_collect )
      if isnothing( covered_indices )
        covering_indices = findall( ref_dentop->ref_dentop⊆dentop, dentop_collect )
        deleteat!( dentop_collect, covering_indices )
        push!( dentop_collect, dentop )
      end # if

      jldopen( shifted_amp_file, "w" ) do shifted_jld_file
        for key ∈ filter( key->key∉("amp_lorentz_list","loop_den_list"), keys(jld_file) )
          shifted_jld_file[key] = jld_file[key]
        end # for key
        shifted_jld_file["amp_lorentz_list"] = to_String( amp_list )
        shifted_jld_file["loop_den_list"] = to_String( den_list )
      end # shifted_jld_file
    else
      den_list = normalize_loop_mom( den_list )
      dentop = DenTop( n_loop, indep_ext_mom, den_list )
      push!( backup_dentop_collect, dentop )
      push!( dentop_collect, dentop )
    end # if

    close( jld_file )
  end # for (index, amp_index)

  unique!( dentop->reduce(*,dentop.den_list), dentop_collect )
  dentop_collect = get_superior_dentop_collect( dentop_collect )
  # @assert all( backup_dentop->any(dentop->backup_dentop⊆dentop,dentop_collect), backup_dentop_collect )
  @info "$(length(dentop_collect)) initial topolgies found."

  @info "Begin to make complete topologies @ $(now())"
  complete_dentop_collect = make_complete_dentop_collect( dentop_collect )
  @assert all( backup_dentop->any(complete_dentop->backup_dentop⊆complete_dentop,complete_dentop_collect), backup_dentop_collect )
  @info "$(length(complete_dentop_collect)) complete topologies found @ $(now())."

  file = open( joinpath( topology_dir, mom_shift_opt ? "shifted_topology.out" : "topology.out" ), "w" )
  line_str = "-"^80
  remaining_indices = (collect∘eachindex)( backup_dentop_collect )
  chosen_indices = Int[]
  for (index, complete_dentop) ∈ enumerate( complete_dentop_collect )
    @assert is_valid_dentop(complete_dentop)
    pos_list = findall( dentop->dentop⊆complete_dentop, backup_dentop_collect )
    isempty( remaining_indices ) && break
    isempty( pos_list ) && continue
    push!( chosen_indices, index )
    numbering = length( chosen_indices )
    setdiff!( remaining_indices, pos_list )

    println()
    println( line_str )
    println( "Complete topology #$(numbering) covers files:" )
    map( println, amp_file_list[ pos_list ] )
    println( line_str )
    map( println, complete_dentop.den_list )

    if mom_shift_opt
      for shifted_amp_file ∈ shifted_amp_file_list[ pos_list ]
        @assert isfile( shifted_amp_file )
        shifted_jld = load( shifted_amp_file )
        if haskey( shifted_jld, "topology_dict" )
          topology_dict = shifted_jld["topology_dict"]
          topology_dict[index] = to_String( complete_dentop.den_list )
          shifted_jld["topology_dict"] = topology_dict
          save( shifted_amp_file, shifted_jld )
        else
          jldopen( shifted_amp_file, "a" ) do shifted_jld
            shifted_jld["topology_dict"] = Dict( index => to_String( complete_dentop.den_list ) )
          end # shifted_jld
        end # if

      end # for shifted_amp_file
    end # if mom_shift_opt

    jldopen( joinpath( topology_dir, mom_shift_opt ? "shifted_topology$(numbering).jld2" : "topology$(numbering).jld2" ), "w" ) do topology_file
      topology_file["n_loop"] = complete_dentop.n_loop
      topology_file["indep_ext_mom"] = to_String( complete_dentop.indep_ext_mom )
      topology_file["den_list"] = to_String( complete_dentop.den_list )

      topology_file["covering_amplitudes"] = mom_shift_opt ? [
        ( amp_file = replace( abspath( amp_file_name ), root_dir => ".." ),
          shifted_amp_file = replace( shifted_amp_file_name, root_dir => ".." ),
          mom_shift = Dict( string(key) => string( mom_shift[key] ) for key ∈ keys(mom_shift) )
        ) for (amp_file_name, shifted_amp_file_name, mom_shift) ∈
          zip( amp_file_list[ pos_list ], shifted_amp_file_list[ pos_list ], mom_shift_collect[ pos_list ] )
      ] : [ replace( abspath( amp_file_name ), root_dir => ".." )
              for amp_file_name ∈ amp_file_list[ pos_list ] ]
    end # topology_file

    cover_amp_str = mom_shift_opt ? 
                      join( [ amp_file_list[pos] *
                              ( isempty(mom_shift_collect[pos]) ? "" :
                                "\n  momenta shift:\n" *
                                join( [ "  ├─ $(key) -> $(value)"
                                        for (key,value) ∈ mom_shift_collect[pos] ],
                                      "\n" ) )
                            for pos ∈ pos_list ], "\n" ) :
                      join( amp_file_list[pos_list], "\n" )
    write( file, """
    $(line_str)
    Complete topology #$(numbering) covers files:
    $(cover_amp_str)
    $(line_str)
    $(join( map(string,complete_dentop.den_list), "\n" ))

    """ )

  end # for (index, complete_dentop)
  close( file )

  @assert isempty(remaining_indices)

  box_message( "Information is in topology.out" )

  return complete_dentop_collect[ chosen_indices ]

end # function construct_den_topology



