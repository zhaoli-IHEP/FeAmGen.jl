

######################################################
const root_dir()::String = (dirname∘dirname∘pathof∘Base.moduleroot)(FeAmGen)
######################################################
const art_dir()::String = Pkg.Artifacts.artifact"FeAmGen"
######################################################



######################################################
"""
    const Dict color_dict

Dict( 1=>`:singlet`, 3=>`:triplet`, -3=>`:triplet`, 8=>`:octet` )
"""
const color_dict = Dict( 1=>:singlet, 3=>:triplet, -3=>:triplet, 8=>:octet )
######################################################
"""
    const Dict spin_dict

Dict( 1=>`:scalar`, -1=>`:ghost`, 2=>`:fermion`, 3=>`:vector` )
"""
const spin_dict = Dict( 1=>:scalar, -1=>:ghost, 2=>:fermion, 3=>:vector )
######################################################
"""
    charge_convert( model_charge::Float64 )::Rational{Int64} 

Convert the float number of charge in the model file into a rational number, especially the quark charge.
"""
charge_convert( model_charge::Float64 )::Rational{Int64} = round(Int64,model_charge*3)//3
######################################################


#########################
"""
    struct Particle

```
  kf::Int64
  name::String
  antiname::String
  spin::Symbol
  color::Symbol
  charge::Rational{Int64}
  mass::Basic
  width::Basic
```
"""
struct Particle
  kf::Int64
  name::String
  antiname::String
  spin::Symbol
  color::Symbol
  charge::Rational{Int64}
  mass::Basic
  width::Basic
end # struct Particle
#########################



#########################
"""
    struct Interaction

```
  name::String
  # UFO convention is outgoing particle, but in our convention is incoming
  link_list::Vector{Particle}
  color_row_list::Vector{Basic} # 1xN
  lorentz_col_list::Vector{Basic} # Mx1
  couplings_matrix::Array{Basic,2} # NxM
  QCD_order::Int64
  QED_order::Int64
  SPC_order::Int64
```
"""
struct Interaction
  name::String
  # UFO convention is outgoing particle, but in our convention is incoming
  link_list::Vector{Particle}
  color_row_list::Vector{Basic} # 1xN
  lorentz_col_list::Vector{Basic} # Mx1
  couplings_matrix::Array{Basic,2} # NxM
  QCD_order::Int64
  QED_order::Int64
  SPC_order::Int64
end # struct Interaction
#########################

#########################
"""
    struct Model

```
  name::String
  unitary_gauge::Bool
  particle_list::Vector{Particle}
  particle_name_dict::Dict{String,Particle}
  particle_kf_dict::Dict{Int64,Particle}
  interaction_list::Vector{Interaction}
  sorted_kf_list_dict::Dict{Vector{Int64},Interaction}
  parameter_dict::Dict{Basic,Basic}
  coupling_dict::Dict{Basic,Basic}
```
"""
struct Model
  name::String
  unitary_gauge::Bool
  particle_list::Vector{Particle}
  particle_name_dict::Dict{String,Particle}
  particle_kf_dict::Dict{Int64,Particle}
  interaction_list::Vector{Interaction}
  sorted_kf_list_dict::Dict{Vector{Int64},Interaction}
  parameter_dict::Dict{Basic,Basic}
  coupling_dict::Dict{Basic,Basic}
end # struct Model
#########################










"""
    is_anticommuting( part::Particle )::Bool

Function "is_anticommuting" is used for the sign in the calculation of fermion/ghost loop.
The problem is according to Paulo Nogueira, we can have our own definition for the sign of majorana loop.
"""
is_anticommuting( part::Particle )::Bool = part.spin in [:fermion, :ghost] 
################################################
"""
    is_neutral( part::Particle )::Bool
"""
is_neutral( part::Particle )::Bool = part.name == part.antiname
"""
    is_not_majorana( part::Particle )::Bool
"""
is_not_majorana( part::Particle )::Bool = part.name != part.antiname
"""
    is_massless( part::Particle )::Bool
"""
is_massless( part::Particle )::Bool = part.mass == 0
"""
    is_massive( part::Particle )::Bool
"""
is_massive( part::Particle )::Bool = part.mass != 0
"""
    is_gluon( part::Particle )::Bool
"""
is_gluon( part::Particle )::Bool = part.kf == 21
"""
    is_photon( part::Particle )::Bool
"""
is_photon( part::Particle )::Bool = part.kf == 22
"""
    is_massless_quark( part::Particle )::Bool
"""
is_massless_quark( part::Particle )::Bool = 1 <= abs(part.kf) <= 5
"""
    is_top_quark( part::Particle )::Bool
"""
is_top_quark( part::Particle )::Bool = abs(part.kf) == 6
"""
    is_ghost( part::Particle )::Bool
"""
is_ghost( part::Particle )::Bool = abs(part.kf) == 9000005
"""
    is_quark( part::Particle )::Bool 
"""
is_quark( part::Particle )::Bool = 1 <= abs(part.kf) <= 6
"""
    is_colorful( part::Particle )::Bool
"""
is_colorful( part::Particle )::Bool = is_gluon(part) || is_quark(part) || is_ghost(part)
"""
    is_not_colorful( part::Particle )::Bool
"""
is_not_colorful( part::Particle )::Bool = !is_colorful(part)
"""
    is_massive_fermion( part::Particle )::Bool
"""
is_massive_fermion( part::Particle )::Bool = part.mass != 0 && part.spin == :fermion
################################################


################################################
"""
    to_string( inter::Interaction )::String

Convert Interaction to string for writing out.
"""
function to_string( inter::Interaction )::String
################################################

  link_name_list = map( p_ -> p_.name, inter.link_list )

  # UFO convention is outgoing particle, but in our convention is incoming
  return inter.name*
    """
      Outgoing: ($( join( link_name_list, "," ) ))
      Color: $( inter.color_row_list )
      Lorentz: $( inter.lorentz_col_list )
      Couplings_matrix: $( inter.couplings_matrix )
      QCD_order: $(inter.QCD_order)
      QED_order: $(inter.QED_order)
      SPC_order: $(inter.SPC_order)
    """
end # function to_string


############################################################################
"""
    find_AC_pair_pos( part_list::Vector{Particle} )::Union{Tuple{Nothing,Nothing},Tuple{Int64,Int64}}

This function is used to find the pair of anti-commuting particles/fields.
This information will be inserted into QGRAF model, so that QGRAF can generate correct sign for anti-commuting particle loop.
For example, we should have extra minus sign for fermion loop and ghost loop.
"""
function find_AC_pair_pos( part_list::Vector{Particle} )::Union{Tuple{Nothing,Nothing},Tuple{Int64,Int64}}
############################################################################

  first_AC_pos = findfirst( p_ -> is_anticommuting(p_), part_list )
  if isnothing(first_AC_pos)
    return nothing, nothing
  end # if

  first_AC_pos = first( first_AC_pos )
  second_AC_pos = findnext( p_ -> is_anticommuting(p_), part_list, first_AC_pos+1 )
  @assert !isnothing(second_AC_pos)
  second_AC_pos = first( second_AC_pos )

  return first_AC_pos, second_AC_pos
end # function find_AC_pair_pos


####################################################################################
"""
    generate_ordered_link_list( part_list::Vector{Particle} )::Vector{Particle}

This function is used to generate order particle list according to the rule of QGRAF.
"""
function generate_ordered_link_list( part_list::Vector{Particle} )::Vector{Particle}
####################################################################################

  first_AC_pos, second_AC_pos = find_AC_pair_pos( part_list )

  clone_link_list = copy( part_list )
  if !isnothing(first_AC_pos) && !isnothing(second_AC_pos)
    first_AC_part = clone_link_list[first_AC_pos]
    @assert first_AC_pos < second_AC_pos
    if first_AC_part.kf > 0 
      clone_link_list[first_AC_pos], clone_link_list[second_AC_pos] = clone_link_list[second_AC_pos], clone_link_list[first_AC_pos]
    end # if
  end # if

  return clone_link_list
end # function generate_ordered_link_list


#############################################################
"""
    to_qgraf_name( name_str::String )::String

Change the particle name into the one that QGRAF can accept.
"""
function to_qgraf_name( name_str::String )::String
#############################################################
  tail_char = name_str[end]
  body_str = name_str[1:end-1]
  if tail_char == '~' 
    return body_str*"bar"
  elseif tail_char == '+'
    return body_str*"plus"
  elseif tail_char == '-'
    return body_str*"minus"
  else
    return name_str
  end # if
end # function to_qgraf_name


################################################
"""
    generate_QGRAF_model( model::Model )::Nothing

Write out the model file that can be used by QGRAF.
"""
function generate_QGRAF_model( model::Model )::Nothing
################################################

  file = open( "model.qgraf", "w" )

  nonnegkf_part_list = filter( p_ -> p_.kf >= 0, model.particle_list )
  sign_str_list = map( p_ -> p_.spin in [:fermion,:ghost] ? "-" : "+", nonnegkf_part_list )
  nonnegkf_part_info_list = map( (p_, s_) -> "[ $(p_.name),$(p_.antiname),$(s_) ]\n", nonnegkf_part_list, sign_str_list )

  write( file, join( nonnegkf_part_info_list, "" ) )
  write( file, "[ QCDct1, QCDct1bar, - ]\n" )
  write( file, "[ QCDct2, QCDct2bar, - ]\n" )


  for inter in model.interaction_list

    ordered_link_list = generate_ordered_link_list( inter.link_list )

    part_name_list = map( p_ -> p_.name, ordered_link_list )
    part_list_str = join( part_name_list, "," )

    n_link = length(inter.link_list)
    if n_link > 2
      base_vertex_str = "[$(part_list_str); epow='$(inter.QED_order)',gspow='$(inter.QCD_order)',spc='$(inter.SPC_order)',qcdct='0']\n"
      write( file, base_vertex_str )
    end # if

    deg_matrix = map( ele_ -> get_degree(ele_,Basic("CTorder")), inter.couplings_matrix )
    max_QCDCT_order = maximum( deg_matrix ) # i.e. maximum degree 
    if max_QCDCT_order == 2 
      QCDct1_vertex_str = "[$(part_list_str), QCDct1bar, QCDct1; epow='$(inter.QED_order)',gspow='$(inter.QCD_order+2)',spc='$(inter.SPC_order)',qcdct='1']\n"
      write( file, QCDct1_vertex_str )
      QCDct2_vertex_str = "[$(part_list_str), QCDct2bar, QCDct2; epow='$(inter.QED_order)',gspow='$(inter.QCD_order+4)',spc='$(inter.SPC_order)',qcdct='2']\n"
      write( file, QCDct2_vertex_str )
    end # if

  end # for iter

  close(file)

  return nothing

end # function generate_QGRAF_model







################################################
"""
    logging_model( model::Model )::Nothing

Write out model file that contains more detail.
"""
function logging_model( model::Model )::Nothing
################################################

  file = open( model.name*".log", "w" )

  write( file, "\n"*"-"^60*"\nAll particles: \n"*"-"^60*"\n" )
  for part in model.particle_list
    write( file, "  ", string(part), "\n" )
  end # for part

  write( file, "\n"*"-"^60*"\nAll interactions: \n"*"-"^60*"\n" )
  for inter in model.interaction_list
    write( file, to_string(inter), "\n" )
  end # for inter

  close(file)

  return nothing

end # function logging_model







###########################################
"""
    vectorized_tensor_product(
      iters...;
      f::Function=identity
    )::Vector{Basic}

Return the vectorized tensor production of several iterable objects.
You could specify what value produced by function `f` applied to the elements.
"""
function vectorized_tensor_product(iters...; f::Function=identity)
###########################################
  return vec( [ f(elem) for elem ∈ Iterators.product(iters...) ] )
end

###########################################
"""
    vectorized_tensor_product_Basic(iters...)
Call `vectorized_tensor_product( iters...; f=pord )` for `Basic`.
"""
vectorized_tensor_product_Basic(iters...) = vectorized_tensor_product(iters...; f=prod)
###########################################
###########################################
"""
    vectorized_tensor_product_Basic(iters...)
Call `vectorized_tensor_product( iters...; f=(elem -> join(elem, ',')) )` for `String`.
"""
vectorized_tensor_product_String(iters...) = vectorized_tensor_product(iters...; f=(elem -> join(elem, ',')))
###########################################

###########################################
function get_diagram_index( file_name::String )::Int
###########################################
  # @assert isfile(amp_file_name)
  file_basename = basename( file_name )

  support_file_type = [ "jld2", "out", "pdf", "tex" ]
  regex_list = map( str -> Regex("[1-9]\\d*\\.$str\$"), support_file_type )

  file_type_index = findfirst( occursin(file_basename), regex_list )
  @assert !isnothing(file_type_index) "Invalid file name: $file_basename"

  str_range = findfirst(regex_list[file_type_index], file_basename)
  str = file_basename[str_range]

  index_str_range = findfirst(r"^[1-9]\d*", str)
  return  parse(Int, str[index_str_range])
end
