###################################
# const definition
const preferred_vac_mom_dict() = Dict{Int,Vector{Vector{Basic}}}(
  1 => [ [ Basic("q1") ] ],
  2 => [ to_Basic( ["q1", "q2", "q1 + q2"] ) ],
  3 => [ to_Basic( ["q1", "q2", "q3", "q1 + q3", "q2 + q3", "q1 + q2 + q3"] ),
         to_Basic( ["q1", "q2", "q3", "q1 + q2", "q1 + q3", "q2 + q3"] ) ],
  4 => [ to_Basic( ["q1", "q2", "q3", "q4",
                    "q1 + q2", "q1 + q3", "q2 + q4",
                    "q1 + q2 + q3", "q1 + q2 + q4",
                    "q1 + q2 + q3 + q4"] ),
         to_Basic( ["q1", "q2", "q3", "q4",
                    "q1 + q2", "q1 + q3", "q2 + q4", "q3 + q4",
                    "q2 + q3 + q4",
                    "q1 + q2 + q3 + q4"] ),
         to_Basic( ["q1", "q2", "q3", "q4",
                    "q1 + q2", "q1 + q3", "q2 + q3", "q2 + q4",
                    "q1 + q2 + q4", "q2 + q3 + q4"] ) ]
) # end preferred_vac_mom_dict
###################################



###########################################
function get_vac_loop_momenta_list(
  ::Val{1} # n_loop
)::Vector{Vector{Basic}}
###########################################

  return [ [ Basic("q1") ] ]

end # function get_vac_loop_momenta_list

###########################################
function get_vac_loop_momenta_list(
  ::Val{2} # n_loop
)::Vector{Vector{Basic}}
###########################################

  return [ to_Basic( ["q1", "q2", "q1 + q2"] )] 

end # function get_vac_loop_momenta_list

###########################################
function get_vac_loop_momenta_list(
  ::Val{3} # n_loop
)::Vector{Vector{Basic}}
###########################################

  return [ to_Basic( ["q1", "q2", "q3", "q1 + q3", "q2 + q3", "q1 + q2 + q3"] ),
           to_Basic( ["q1", "q2", "q3", "q1 + q2", "q1 + q3", "q2 + q3"] ) ]

end # function get_vac_loop_momenta_list

###########################################
function get_vac_loop_momenta_list(
  ::Val{4} # n_loop
)::Vector{Vector{Basic}}
###########################################

  return [ to_Basic( ["q1", "q2", "q3", "q4",
                     "q1 + q2", "q1 + q3", "q2 + q4",
                     "q1 + q2 + q3", "q1 + q2 + q4",
                     "q1 + q2 + q3 + q4"] ),
          to_Basic( ["q1", "q2", "q3", "q4",
                     "q1 + q2", "q1 + q3", "q2 + q4", "q3 + q4",
                     "q2 + q3 + q4",
                     "q1 + q2 + q3 + q4"] ),
          to_Basic( ["q1", "q2", "q3", "q4",
                     "q1 + q2", "q1 + q3", "q2 + q3", "q2 + q4",
                     "q1 + q2 + q4", "q2 + q3 + q4"] ) ]
           
end # function get_vac_loop_momenta_list




##############################################
# den_list: [ Den(...), Den(...), ... ]
function is_matched_preferred_vac_mom_list(
    den_list::Vector{Basic}, 
)::Bool
##############################################

  mom_list = map( x->(first∘get_args)(x), den_list )
  symbol_list = free_symbols(mom_list)
  qi_list = filter( is_loop_mom, symbol_list )

  ext_mom_list = setdiff( symbol_list, qi_list )
  vanish_map = Dict{Basic,Basic}( ext_mom_list .=> zero(Basic) )
  input_vac_mom_list = (unique∘map)( x->subs(x,vanish_map), mom_list )

  n_loop = length(qi_list)
  vac_mom_list_collection = preferred_vac_mom_dict()[n_loop]
  for vac_mom_list in vac_mom_list_collection
    if (isempty∘setdiff)( input_vac_mom_list, vac_mom_list )
      return true
    end # if
  end # for vac_mom_list

  return false

end # function is_matched_preferred_vac_mom_list
















###################################
# This function is used for the case before canonicalization, 
#   so we only check the coefficient sign of the leading qi.
function normalize_loop_mom_single( 
    loop_mom::Basic
)::Basic
###################################

  leading_qi = (first∘get_loop_momenta)( [ loop_mom ] )
  # qi_list = free_symbols(loop_mom)
  # filter!(sym -> (first∘string)(sym) == 'q', qi_list)
  # sort!(qi_list, by=q->parse(Int, string(q)[2:end]))
  # leading_qi = first(qi_list)

  loop_mom = expand(loop_mom)
  if SymEngine.coeff(loop_mom,leading_qi) < 0
    return expand(-loop_mom)
  end # if

  return loop_mom

end # function normalize_loop_mom_single


###################################
function normalize_loop_mom( 
    loop_den_list::Vector{Basic} 
)::Vector{Basic}
###################################

  @funs Den

  qi_list = get_loop_momenta(loop_den_list)
  # qi_list = free_symbols(loop_den_list)
  # filter!(sym -> (first∘string)(sym) == 'q', qi_list)
  # sort!(qi_list, by=q->parse(Int, string(q)[2:end]))

  new_loop_den_list = zeros( Basic, length(loop_den_list) )
  for index in eachindex(loop_den_list)
    mom, mass, width = get_args(loop_den_list[index])
    mom = expand(mom)
    coeff_list = SymEngine.coeff.( mom, qi_list )
    first_nonzero_index = findfirst( !iszero, coeff_list )
    @assert !isnothing(first_nonzero_index)
    if coeff_list[first_nonzero_index] < 0
      mom = expand(-mom)
    end # if
    new_loop_den_list[index] = Den( mom, mass, width )
  end # for index

  return new_loop_den_list

end # function normalize_loop_mom






#########################################################
# Created by Quan-feng Wu
# Feb. 16 2023
function gen_loop_mom_canon_map(
    mom_list::Vector{Basic},
    mom_list_collect::Vector{Vector{Basic}}=Vector{Basic}[]
)::Dict{Basic,Basic}
#########################################################

  # get information
  q_list = get_loop_momenta( mom_list )
  k_list = get_ext_momenta( mom_list )
  n_loop = isempty(q_list) ? 0 : (get_loop_index∘last)( q_list )
  @assert q_list == [Basic("q$ii") for ii ∈ 1:n_loop]
  q_null_dict = Dict( q_list .=> zero(Basic) )
  k_null_dict = Dict( k_list .=> zero(Basic) )
  preferred_flag = n_loop ∈ keys(preferred_vac_mom_dict())

  # default repl_rule
  chosen_repl_rule = Dict{Basic,Basic}(q_list .=> q_list)
  chosen_repl_order = Inf

  # find all branches of the loop momenta.
  mom_q_coeff_mat = coefficient_matrix( mom_list, q_list )
  @assert (isempty∘setdiff)( unique(mom_q_coeff_mat), Basic[-1,0,1] )
  mom_q_coeff_list = filter( !iszero, (collect∘eachrow)( mom_q_coeff_mat ) )
  mom_q_coeff_list = map( coeff_list->coeff_list./coeff_list[findfirst(!iszero,coeff_list)],
                            mom_q_coeff_list )
  unique!( mom_q_coeff_list )
  sort!( mom_q_coeff_list; 
          by=coeff_list->begin
            num_q = (length∘findall)(!iszero,coeff_list)
            new_coeff_list = map( coeff->subs(coeff,Dict(Basic(-1)=>2,Basic(0)=>3)), coeff_list )
            num_q, new_coeff_list
          end )
  vac_mom_list = [ sum( mom_q_coeff .* q_list ) for mom_q_coeff ∈ mom_q_coeff_list ]
  branch_indices_list = Vector{Int}[]
  for mom_q_coeff ∈ mom_q_coeff_list
    branch_indices = findall( row->row==mom_q_coeff||row==-mom_q_coeff,
                                eachrow(mom_q_coeff_mat) )
    push!( branch_indices_list, branch_indices )
  end # for mom_q_coeff

  # cover flag
  chosen_cover_flag = false
  has_cover_flag = false

  for selected_branch_indices_list ∈ permutations( branch_indices_list, length(q_list) ),
      selected_mom_indices ∈ Iterators.product( selected_branch_indices_list... )
    for sign_list ∈ Iterators.product([(1, -1) for _ in q_list]...)

      # check invertible q coeff matrix
      selected_vac_mom_list = map( mom->subs(mom,k_null_dict), mom_list[collect(selected_mom_indices)] )
      selected_coeff_mat = coefficient_matrix( selected_vac_mom_list, q_list )
      (iszero∘expand∘get_det)( selected_coeff_mat ) && break
      inv_selected_coeff_mat = inv( selected_coeff_mat )

      # check q signs are all positive or all negative
      vac_repl_rule = Dict( q_list .=> inv_selected_coeff_mat * (q_list .* sign_list) )
      new_vac_mom_list = map( vac_mom->(expand∘subs)(vac_mom,vac_repl_rule), vac_mom_list )
      map!( normalize_loop_mom_single, new_vac_mom_list, new_vac_mom_list )
      new_coeff_mat = coefficient_matrix( new_vac_mom_list, q_list )
      !all( ≥(0), new_coeff_mat ) && continue

      # check preferred vacumm momenta configuration
      preferred_flag &&
        !any( preferred_mom_list->new_vac_mom_list⊆preferred_mom_list,
                preferred_vac_mom_dict()[n_loop] ) && break

      # construct the repl rule
      selected_ext_mom_list = map( mom->subs(mom,q_null_dict), mom_list[collect(selected_mom_indices)] )
      this_repl_rule = Dict( q_list .=>
                              expand.( inv_selected_coeff_mat * (sign_list .* q_list .- selected_ext_mom_list) ) )
      this_mom_list = map( mom->(expand∘subs)(mom,this_repl_rule), mom_list )
      map!( normalize_loop_mom_single, this_mom_list, this_mom_list )
      (!isempty∘setdiff)( 
        (unique∘map)(abs, coefficient_matrix( this_mom_list, k_list )), Basic[0,1] ) && continue

      # check cover
      this_cover_flag = false
      if any( mom_list->this_mom_list⊆mom_list||mom_list⊆this_mom_list, mom_list_collect )
        has_cover_flag = true
        this_cover_flag = true
      elseif has_cover_flag
        continue
      end # if
      this_repl_order = get_sort_order( this_mom_list )

      # check order
      if this_cover_flag && !chosen_cover_flag
        chosen_repl_order = this_repl_order
        chosen_repl_rule = this_repl_rule
        chosen_cover_flag = true
      elseif this_repl_order < chosen_repl_order
        chosen_repl_order = this_repl_order
        chosen_repl_rule = this_repl_rule
      # elseif this_repl_order == chosen_repl_order
        # @warn "Two equivalent replacement rules are found."
      end # if
    end # for sign_list
  end # for selected_branch_indices_list, selected_mom_indices

  return chosen_repl_rule

end # function gen_loop_mom_canon_map






#########################################################
# Created by Quan-feng Wu 
# Mar. 26 2023
# 
# A simple trial for sorting the replace rules generated in `gen_loop_mom_canon_map`.
function get_sort_order( 
    mom_list::Vector{Basic} 
)::BigInt
#########################################################

  tmp_mom_list = unique( abs, mom_list )
  map!( normalize_loop_mom_single, tmp_mom_list, tmp_mom_list )
  q_list = get_loop_momenta( tmp_mom_list )
  k_list = get_ext_momenta( tmp_mom_list )
  qk_list = vcat(q_list, k_list)

  order = one(BigInt)

  for (mom_index, mom) ∈ enumerate(tmp_mom_list)
    mom_qk_coeff_mat = coefficient_matrix( [mom], qk_list )

    qk_str = (join∘map)( coeff->coeff==-1 ? "2" : string(coeff), mom_qk_coeff_mat )

    order *= parse(BigInt, reverse(qk_str), base=3)
  end # for (mom_index, mom)

  return order
end # function get_sort_order






#####################################################
function canonicalize_amp( 
    loop_den_list::Vector{Basic},  
    amp_lorentz_list::Vector{Basic}
    # mom_list_collect::Vector{Vector{Basic}}=Vector{Basic}[]
)::Tuple{Vector{Basic},Vector{Basic}}
######################################################

  # n_loop = get_n_loop( loop_den_list )
  q_list = get_loop_momenta( loop_den_list )
  # k_list = get_ext_momenta( loop_den_list )
  n_loop = isempty(q_list) ? 0 : (get_loop_index∘last)( q_list )

  n_loop == 0 && return loop_den_list, amp_lorentz_list 

  mom_list = map( first∘get_args, loop_den_list )
  canon_map = gen_loop_mom_canon_map( mom_list )
  # canon_map = gen_loop_mom_canon_map( mom_list, mom_list_collect )

  new_loop_den_list = map( den->subs(den,canon_map), loop_den_list )
  new_loop_den_list = normalize_loop_mom( new_loop_den_list )
  new_amp_lorentz_list = map( amp->subs(amp,canon_map), amp_lorentz_list )

  new_mom_list = map( first∘get_args, new_loop_den_list )

  # CHECK begin
  # qi_list = Basic[ Basic("q$ii") for ii in 1:n_loop ]
  unique_coeff_list = union( coefficient_matrix( new_mom_list, q_list ) )
  @assert all( ≥(0), unique_coeff_list ) "$new_mom_list"
  # CHECK end

  return new_loop_den_list, new_amp_lorentz_list 

end # function canonicalize_amp


