using AbstractAlgebra
using JLD2
using SymEngine
using FeynUtils


##################################
function get_coeff_row_mat(
    one_den::Basic, 
    sp_dict::Dict{Basic,Basic}, 
    sp_var_list::Vector{Basic} 
)::Matrix{Rational}
##################################

  one_mom = (expand∘first∘get_args)(one_den)
  mom2 = subs( make_SP(one_mom,one_mom), sp_dict )

  n_sp = length(sp_var_list)
  coeff_row_mat = zeros( Rational, 1, n_sp )
  for sp_index in 1:n_sp
    sp_coeff = SymEngine.coeff( mom2, sp_var_list[sp_index] )
    coeff_row_mat[1,sp_index] = Rational( parse( Int64, string(sp_coeff) ) )
  end # for sp_index

  return coeff_row_mat

end # function get_coeff_row_mat


##################################
function find_indep_den_list(
    den_list::Vector{Basic}, 
    den_universe::Vector{Basic} 
)::Vector{Basic}
##################################

  mom_list = (free_symbols∘map)( (first∘get_args), den_universe )
  qi_list = filter( x->(first∘string)(x)=='q', mom_list )
  indep_mom_list = setdiff( mom_list, qi_list )
  n_loop = length(qi_list)

  sp_index = 1::Int64
  sp_var_list = Vector{Basic}()
  sp_dict = Dict{Basic,Basic}()

  for i1 in 1:n_loop, i2 in i1:n_loop
    qi = qi_list[i1] 
    qj = qi_list[i2]
    sp_var = Basic("sp$(sp_index)")
    push!( sp_var_list, sp_var )
    push!( sp_dict, make_SP(qi,qj) => sp_var )
    sp_index += 1
  end # for qi_index, qj_index

  for one_mom in indep_mom_list, qi_index in 1:n_loop
    qi = qi_list[qi_index]
    sp_var = Basic( "sp$(sp_index)" )
    push!( sp_var_list, sp_var )
    push!( sp_dict, make_SP(one_mom,qi) => sp_var )
    sp_index += 1
  end # for one_mom

  n_sp = length(sp_var_list)
  coeff_mat = zeros( Rational, 0, n_sp )
  for one_den in den_list
    new_row_mat = get_coeff_row_mat( one_den, sp_dict, sp_var_list )

    new_coeff_mat = vcat( coeff_mat, new_row_mat )
    @assert rank(new_coeff_mat) > rank(coeff_mat)
    coeff_mat = new_coeff_mat
  end # for one_den 
  rank_coeff_mat = rank(coeff_mat)

  indep_den_list = Vector{Basic}()
  for one_den in den_universe
    trial_row_mat = get_coeff_row_mat( one_den, sp_dict, sp_var_list )

    if (rank∘vcat)( coeff_mat, trial_row_mat ) > rank_coeff_mat
      push!( indep_den_list, one_den )
    end # if
  end # for one_den

  return indep_den_list

end # function find_indep_den_list





###############################################
function find_first_parent(
    dentop::Basic, # DT(...)
    dentop_collect::Vector{Basic} # [ DT(...), DT(...), ... ] 
)::Union{Int64,Nothing}
###############################################

  # Search for the parent of den_list
  pos = findfirst( x -> get_args(dentop) ⊆ get_args(x), dentop_collect )
  return pos

end # function find_first_parent


###############################################
function find_all_children(
    dentop::Basic, # DT(...)
    dentop_collect::Vector{Basic} # [ DT(...), DT(...), ... ]
)::Vector{Int64}
###############################################

  # Search for the parent of den_list
  pos_list = findall( x -> get_args(x) ⊆ get_args(dentop), dentop_collect )
  return pos_list

end # function find_all_children


#############################################
function get_superior_dentop_chain(
    dentop_chain::Basic # DT*DT*...
)::Basic 
#############################################

  dentop_collect = get_mul_vector(dentop_chain)

  new_dentop_collect = Vector{Basic}()
  for dentop in dentop_collect

    pos = find_first_parent( dentop, new_dentop_collect )
    if pos != nothing 
      # Found parent and check next.
      continue
    end # if

    # Remove all children and then insert as parent.
    pos_list = find_all_children( dentop, new_dentop_collect )
    for pos in pos_list
      new_dentop_collect[pos] = one(Basic)
    end # for pos

    push!( new_dentop_collect, dentop )
  end # for one_file 
  #new_dentop_collect = filter( !isone, new_dentop_collect )

  return reduce( *, new_dentop_collect )

end # function get_superior_dentop_chain




#####################################################
function get_top_cover_combination(
    dentop_chain::Basic, # DT*DT*...
    den_universe::Vector{Basic} 
)::Tuple{Int64,Basic} # top_cover, DT*DT*...+DT*DT*...+....
#####################################################

  @funs DT

  # input: dentop_chain, den_universe
  # Iterate over all dentop.
  max_cover_dict = Dict{Basic,Int64}()

  dentop_collect = get_mul_vector(dentop_chain)
  for dentop in dentop_collect
    den_list = get_args(dentop)
    # Find set of Den(...) in loop_den_market that is independent of den_list.
    indep_den_list = find_indep_den_list( den_list, den_universe )

    # Find cover size for each extra_den in indep_den_list.
    cover_dict = Dict{Basic,Int64}()
    for extra_den in indep_den_list
      trial_dentop = DT( den_list..., extra_den )
      pos_list = find_all_children( trial_dentop, dentop_collect )
      cover_size = length(pos_list)-1 # den_list is always child 
      if iszero(cover_size) 
        continue
      end # if
      push!( cover_dict, extra_den => cover_size )
    end # for extra_den

    if isempty(cover_dict)
      # This den_list cannot make any cover via only one extension.
      continue
    end # if

    max_cover = max( values(cover_dict)... )
    extra_den_list = (collect∘keys∘filter)( x->last(x)==max_cover, cover_dict )

    # possibilities are summed up DT+DT+... for max_cover
    new_DT_sum = (sum∘map)( x -> DT(den_list...,x), extra_den_list )
    push!( max_cover_dict, new_DT_sum => max_cover )
  end # for den_list

  max_cover_list = (collect∘values)(max_cover_dict) 
  top_cover = isempty(max_cover_list) ? 0 : max( max_cover_list... )

  top_cover_combination = one(Basic)
  for (new_DT_sum,max_cover) in max_cover_dict
    if max_cover != top_cover
      continue
    end # if
    top_cover_combination *= new_DT_sum
  end # for one_pair

  return top_cover, expand(top_cover_combination)

end # function get_top_cover_combination





##############################################################
function get_new_multiverse(
    dentop_chain_multiverse::Vector{Basic}, # [ DT*DT*..., DT*DT*..., ... ]
    den_universe::Vector{Basic} # [ Den(...), Den(...), ... ]
)::Tuple{ Vector{Basic}, Vector{Basic} }
##############################################################

  new_dentop_chain_multiverse = Vector{Basic}()
  complete_dentop_chain_multiverse = Vector{Basic}()
  for one_dentop_chain in dentop_chain_multiverse

    # Find the max cover by combining 
    #   one dentop from dentop_chain 
    #   and one den from den_universe. 
    # top_cover, DT*DT*...+DT*DT*...+...
cost_time = @elapsed begin
    top_cover, top_cover_combination = get_top_cover_combination( one_dentop_chain, den_universe )
end # cost_time
println( "cost $(cost_time) sec" )

    if isone(top_cover_combination)
      push!( complete_dentop_chain_multiverse, one_dentop_chain )
      continue
    end # if

    # implement combination
    new_split_multiverse = get_add_vector_expand(one_dentop_chain*top_cover_combination)
cost_time = @elapsed begin
    for index in 1:length(new_split_multiverse)
      new_split_multiverse[index] = get_superior_dentop_chain( new_split_multiverse[index] )
    end # for one_multiverse
end # cost_time
println( "update multiverse cost $(cost_time) sec" )
    append!( new_dentop_chain_multiverse, new_split_multiverse )

  end # for one_dentop_chain

  return new_dentop_chain_multiverse, complete_dentop_chain_multiverse

end # function get_new_multiverse






###########################
function main()::Nothing
###########################

  @funs DT

  dir = "b_g_TO_Wminus_t_2Loop_amplitudes"
  root, dirs, files = (first∘collect∘walkdir)(dir)
  file_list = filter( s->endswith(s,".jld2"), files )

  dentop_collect = Vector{Basic}() # [ DT(...), DT(...), ... ]
  for one_file in file_list
    file = jldopen( "$(dir)/$(one_file)", "r" )
    den_list = (to_Basic∘read)( file, "loop_den_list" )
    close( file )

    push!( dentop_collect, DT(den_list...) )
  end # for one_file
  unique!( dentop_collect )
  dentop_chain = reduce( *, dentop_collect )
@show length(dentop_collect)
error("DEBUG")

  # DT(...)*DT(...)*...
  dentop_chain = get_superior_dentop_chain( dentop_chain )
  @show length(file_list) (length∘get_mul_vector)(dentop_chain) 

  # Create den_universe.
  den_universe = union( map( x->map(get_args,get_mul_vector(x)), dentop_chain )... )
  @show den_universe length(den_universe)

  dentop_chain_multiverse = Basic[ dentop_chain ] 
  complete_dentop_chain_multiverse = Vector{Basic}() 

  while !isempty( dentop_chain_multiverse )
    println( "Try to generate new branches" )
    new_dentop_chain_multiverse, new_complete_dentop_chain_multiverse = 
        get_new_multiverse( dentop_chain_multiverse, den_universe )
    println( "Got it" )
    complete_dentop_chain_multiverse = vcat( complete_dentop_chain_multiverse, 
                                           new_complete_dentop_chain_multiverse ) 

@show length(dentop_chain_multiverse)

    dentop_chain_multiverse = new_dentop_chain_multiverse

@show length(dentop_chain_multiverse)

    min_length = min( map( (length∘get_mul_vector), dentop_chain_multiverse )... )
@show min_length
    dentop_chain_multiverse = filter( x->(length∘get_mul_vector)(x)==min_length, dentop_chain_multiverse )

    n_branch = length(dentop_chain_multiverse)
@show n_branch
    for index in 1:n_branch
#      println( "Branch #$(index)" )
      new_dentop_chain = dentop_chain_multiverse[index]
##    println( (length∘get_mul_vector)(new_dentop_chain) ) 
    end # for index
  end # while


  return nothing

end # function main

#########
main()
#########

