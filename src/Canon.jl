###################################
# const definition
const preferred_vac_mom_Dict() = Dict{Int,Vector{Vector{Basic}}}(
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
) # end preferred_vac_mom_Dict
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
    loop_den_mom_list_collect::Vector{Vector{Basic}}=Vector{Basic}[]
)::Dict{Basic,Basic}
#########################################################

  q_list = get_loop_momenta( mom_list )
  k_list = get_ext_momenta( mom_list )
  n_loop = isempty(q_list) ? 0 : (get_loop_index∘last)( q_list )
  @assert q_list == [Basic("q$ii") for ii ∈ 1:n_loop]

  # Find all branches of the loop momenta.
  mom_q_coeff_mat = coefficient_matrix( mom_list, q_list )
  mom_q_coeff_list = unique( row->sort([row,-row]), eachrow(mom_q_coeff_mat) )
  filter( !iszero, mom_q_coeff_list )
  branch_indices_list = Vector{Int}[]
  for mom_q_coeff ∈ mom_q_coeff_list
    branch_indices = findall( row->row==mom_q_coeff||row==-mom_q_coeff,
                                eachrow(mom_q_coeff_mat) )
    push!( branch_indices_list, branch_indices )
  end # for mom_q_coeff

  vac_mom_list = map( mom_q_coeff->transpose(mom_q_coeff)*q_list, mom_q_coeff_list )

  preferred_flag = n_loop ∈ keys(preferred_vac_mom_Dict())

  chosen_repl_order = Inf
  chosen_repl_rule = Dict{Basic,Basic}()
  chosen_cover_flag = false
  has_cover_flag = false

  for selected_vac_mom_indices ∈ permutations( eachindex(vac_mom_list), length(q_list) )
    for sign_list ∈ Iterators.product([(1, -1) for _ in q_list]...)
      last(sign_list) == -1 && break

      selected_vac_mom_list = map( expand, sign_list .* vac_mom_list[selected_vac_mom_indices] )
      selected_coeff_mat = coefficient_matrix( selected_vac_mom_list, q_list )
      (iszero∘expand∘get_det)( selected_coeff_mat ) && break

      vac_repl_rule = Dict( q_list .=> inv( selected_coeff_mat ) * q_list )
      new_vac_mom_list = map( vac_mom->(expand∘subs)(vac_mom,vac_repl_rule), vac_mom_list )
      @assert new_vac_mom_list[selected_vac_mom_indices] == sign_list .* q_list

      map!( normalize_loop_mom_single, new_vac_mom_list, new_vac_mom_list )
      new_coeff_mat = coefficient_matrix( new_vac_mom_list, q_list )
      !all( ≥(0), new_coeff_mat ) && continue

      preferred_flag &&
        !any( preferred_mom_list->new_vac_mom_list⊆preferred_mom_list,
                preferred_vac_mom_Dict()[n_loop] ) && break

      for mom_indices ∈ Iterators.product( branch_indices_list[selected_vac_mom_indices]... )
        tmp_k_coeff_mat = coefficient_matrix( mom_list[collect(mom_indices)], k_list )
        this_repl_rule = Dict{Basic,Basic}()
        new_repl = Dict{Basic,Basic}( q_list[q_index] => q_list[q_index] -
                                        transpose(tmp_k_coeff_mat[ q_index, : ]) * k_list
                                          for q_index ∈ eachindex(q_list) )
        for key ∈ keys(vac_repl_rule)
          this_repl_rule[key] = (expand∘subs)( vac_repl_rule[key], new_repl )
        end # for key

        tmp_mom_list = map( mom->(expand∘subs)(mom,this_repl_rule), mom_list )
        (!isempty∘setdiff)( 
          (union∘map)(abs, coefficient_matrix( tmp_mom_list, k_list )), Basic[0,1] ) && continue
        tmp_cover_flag = false

        if any( mom_list->tmp_mom_list⊆mom_list||mom_list⊆tmp_mom_list, loop_den_mom_list_collect )
          has_cover_flag = true
          tmp_cover_flag = true
        elseif has_cover_flag
          continue
        end # if
        tmp_repl_order = get_sort_order( tmp_mom_list )

        if tmp_cover_flag && !chosen_cover_flag
          chosen_repl_order = tmp_repl_order
          chosen_repl_rule = this_repl_rule
          chosen_cover_flag = true
        elseif tmp_repl_order < chosen_repl_order
          chosen_repl_order = tmp_repl_order
          chosen_repl_rule = this_repl_rule
          # @show chosen_repl_order
        elseif tmp_repl_order == chosen_repl_order
          #@warn "Two equivalent replacement rules are found."
          # _, index = findmin( string, [tmp_repl_order, chosen_repl_order] )
          # index == 1 && (chosen_repl_rule = this_repl_rule)
        end # if
      end # for mom_indices
    end # for sign_list
  end # for selected_mom_indices

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
  # q_index_list = map( get_loop_index, q_list )
  k_list = get_ext_momenta( tmp_mom_list )
  # k_index_list = map( get_ext_index, k_list )
  qk_list = vcat(q_list, k_list)

  # order_vec = Vector{String}(undef, length(tmp_mom_list))
  order = one(BigInt)

  for (mom_index, mom) ∈ enumerate(tmp_mom_list)
    mom_qk_coeff_mat = coefficient_matrix( [mom], qk_list )
    # mom_k_coeff_mat = coefficient_matrix( [mom], k_list )

    qk_str = (join∘map)( coeff->coeff==-1 ? "2" : string(coeff), mom_qk_coeff_mat )
    # k_str = (join∘map)( coeff->coeff==-1 ? "2" : string(coeff), mom_k_coeff_mat )
    # @show mom, qk_str

    order *= parse(BigInt, reverse(qk_str), base=3)
  end # for (mom_index, mom)

  # sort!( order_vec; by=string )

  return order
end # function get_sort_order






#####################################################
function canonicalize_amp( 
    loop_den_list::Vector{Basic},  
    amp_lorentz_list::Vector{Basic},
    loop_den_mom_list_collect::Vector{Vector{Basic}}=Vector{Basic}[]
)::Tuple{Vector{Basic},Vector{Basic}}
######################################################

  # n_loop = get_n_loop( loop_den_list )
  q_list = get_loop_momenta( loop_den_list )
  # k_list = get_ext_momenta( loop_den_list )
  n_loop = isempty(q_list) ? 0 : (get_loop_index∘last)( q_list )

  n_loop == 0 && return loop_den_list, amp_lorentz_list 

  mom_list = map( first∘get_args, loop_den_list )
  canon_map = gen_loop_mom_canon_map( mom_list, loop_den_mom_list_collect ) 

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


