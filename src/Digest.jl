
###########################################################
"""
    check_has_yt( 
        coupling_dict::Dict{Any,Any} 
    )::Bool

Check if there is "yt" in the values of couplings in one "Vertex".
"""
function check_has_yt( 
    coupling_dict::Dict{Any,Any} 
)::Bool
###########################################################

  value_list_str = join( map( v_ -> v_.value, values(coupling_dict) ), "," )
  has_yt = occursin( "yt", value_list_str )

  return has_yt

end # function check_has_yt



#########################################################################
"""
    extract_QCD_QED_order( 
        coupling_dict::Dict{Any,Any} 
    )::Tuple{Int64,Int64,Int64}

Extract the QCD, QED and special orders from the python model file 
  for the "couplings" property in one "Vertex".
"""
function extract_QCD_QED_order( 
    coupling_dict::Dict{Any,Any} 
)::Tuple{Int64,Int64,Int64}
##########################################################################

  QCD_order_set = Set{Int64}( map( v_ -> haskey(v_.order,"QCD") ? v_.order["QCD"] : 0, values(coupling_dict) ) )
  QED_order_set = Set{Int64}( map( v_ -> haskey(v_.order,"QED") ? v_.order["QED"] : 0, values(coupling_dict) ) )
  HIG_order_set = Set{Int64}( map( v_ -> haskey(v_.order,"HIG") ? v_.order["HIG"] : 0, values(coupling_dict) ) )
  SPC_order_set = Set{Int64}( map( v_ -> haskey(v_.order,"SPC") ? v_.order["SPC"] : 0, values(coupling_dict) ) )

  @assert length(QCD_order_set) == 1
  QCD_order = first(QCD_order_set)
  @assert length(QED_order_set) == 1
  QED_order = first(QED_order_set)
  @assert length(HIG_order_set) == 1
  HIG_order = first(HIG_order_set)
  @assert length(SPC_order_set) == 1
  SPC_order = first(SPC_order_set)

  return QCD_order, max(QED_order,HIG_order), SPC_order

end # function extract_QCD_QED_order


###########################################################################
"""
    calculate_CTcoeff( 
        link_part_list::Vector{Particle}, 
        has_yt::Bool, 
        QCD_order::Int64 
    )::Basic

Generate the expression for the QCD counter-term 
  coefficient up-to second order for this vertex.
"""
function calculate_CTcoeff( 
    link_part_list::Vector{Particle}, 
    has_yt::Bool, 
    QCD_order::Int64 
)::Basic
############################################################################

  local CTcoeff = Basic(" (1+dZgx1*CTorder+dZgx2*CTorder^2)^$QCD_order ")
  if has_yt == true 
    CTcoeff *= Basic(" 1+dZmtx1*CTorder+dZmtx2*CTorder^2 ")
  end # if

  for link_part in link_part_list
    if is_gluon(link_part) 
      CTcoeff *= Basic(" sqrt( 1+dZ3x1*CTorder+dZ3x2*CTorder^2 ) ")
    elseif is_massless_quark(link_part)
      CTcoeff *= Basic(" sqrt( 1+dZ2x1*CTorder+dZ2x2*CTorder^2 ) ")
    elseif is_top_quark(link_part)
      CTcoeff *= Basic(" sqrt( 1+dZ2tx1*CTorder+dZ2tx2*CTorder^2 ) ")
    elseif is_ghost(link_part)
      CTcoeff *= Basic(" sqrt( 1+dtZ3x1*CTorder+dtZ3x2*CTorder^2 ) ")
    else
      @assert is_not_colorful(link_part)
    end # if
  end # for link_part 

  CTcoeff = series( CTcoeff, SymEngine.symbols("CTorder"), 0, 2 )

  return CTcoeff

end # function calculate_CTcoeff



#######################################################################################
"""
    extract_couplings_matrix( 
        couplings_dict::Dict{Any,Any}, 
        color_dim::Int64, 
        lorentz_dim::Int64, 
        CTcoeff::Basic 
    )::Array{Basic,2}

Translate the coupling Dict into a two-dimensional matrix, which can be used later.
"""
function extract_couplings_matrix( 
    couplings_dict::Dict{Any,Any}, 
    color_dim::Int64, 
    lorentz_dim::Int64, 
    CTcoeff::Basic 
)::Array{Basic,2}
#########################################################################################

  local couplings_matrix = zeros(Basic,color_dim,lorentz_dim)
  for ele_coup in couplings_dict
    key = ele_coup[1]
    value_str = lowercase( ele_coup[2].name )
    # Since here key starts from (0,0) and in julia it should be (1,1)
    couplings_matrix[key.+1...] = Basic( value_str )*CTcoeff
  end # for ele_coup

  return couplings_matrix

end # function extract_couplings_matrix









#########################################################################
"""
    generate_color_lorentz_couplings( 
        part::Particle 
    )::Tuple{Array{Basic,1},Array{Basic,1},Array{Basic,2}}

Manually generate the color and lorentz couplings for selfenergy line, 
  which is indicated by `part::Particle`.
"""
function generate_color_lorentz_couplings( 
    part::Particle 
)::Tuple{Array{Basic,1},Array{Basic,1},Array{Basic,2}}
###########################################################################

  color_row_list = [ Basic("Identity(1,2)") ]

  mass_str = string(part.mass)

  if is_massless_quark(part) 
    lorentz_col_list = [ Basic("Gamma(-1,1,2)*P(-1,1)") ]
    couplings_matrix = zeros(Basic,1,1)
    couplings_matrix[1,1] = Basic("(-I)*( dZ2x1*CTorder+dZ2x2*CTorder^2 )")
  elseif is_top_quark(part)
    lorentz_col_list = [ Basic("Gamma(-1,1,2)*P(-1,1)"), Basic("Identity(2,1)") ]
    couplings_matrix = zeros(Basic,1,2)
    couplings_matrix[1,1] = Basic("I*( dZ2tx1*CTorder+dZ2tx2*CTorder^2 )")
    couplings_matrix[1,2] = Basic("(-I)*$(mass_str)*( (dZ2tx1+dZmtx1)*CTorder+(dZ2tx1*dZmtx1+dZ2tx2+dZmtx2)*CTorder^2 )")
  elseif is_gluon(part)
    lorentz_col_list = [ Basic("P(1,1)*P(2,1)"), Basic("Metric(1,2)*P(-1,1)*P(-1,1)") ]
    couplings_matrix = zeros(Basic,1,2)
    couplings_matrix[1,1] = Basic("I*( dZ3x1*CTorder+dZ3x2*CTorder^2 )")
    couplings_matrix[1,2] = Basic("(-I)*( dZ3x1*CTorder+dZ3x2*CTorder^2 )")
  elseif is_ghost(part)
    lorentz_col_list = [ Basic("P(-1,1)*P(-1,1)") ]
    couplings_matrix = zeros(Basic,1,1)
    couplings_matrix[1,1] = Basic("I*( dtZ3x1*CTorder+dtZ3x2*CTorder^2 )")
  else
    @assert false "Exception for particle: "*part.name
  end # if

  return color_row_list, lorentz_col_list, couplings_matrix

end # function generate_color_lorentz_couplings










#############################################################
"""
    readin_model( input::Dict{Any,Any} )::Model

Read-in the model detail from python model file.
In this seed program, we only need particle list.
The directory of model files are supposed in ".".
"""
function readin_model(
  input::Dict{Any,Any};
  model_paths::Vector{String}=[pwd()]
)::Model
#############################################################

  model_name = input["model_name"]

  # append the path that python can find the model files
  sys = pyimport( "sys" )
  model_path_index = findfirst( isfile, joinpath.( model_paths, model_name, "object_library.py" ) )
  if isnothing(model_path_index)
    if (!isfile∘joinpath)( art_dir(), "Models", model_name, "object_library.py" )
      error("Model $(model_name) not found @ $model_paths and \"$(art_dir())/Models\".")
    end
    push!( sys."path", "$(art_dir())/Models" )
  else
    push!( sys."path", model_paths[model_path_index] )
  end

  # For example sm.object_library include the basic structure of this model
  py_model = pyimport( "$(model_name).object_library" )

  particle_list = Vector{Particle}()
  particle_name_dict = Dict{String,Particle}()
  particle_kf_dict = Dict{Int64,Particle}()
  for part in py_model.all_particles
    new_part = Particle( part.pdg_code, to_qgraf_name(part.name), to_qgraf_name(part.antiname),
                         spin_dict[part.spin], color_dict[part.color], charge_convert(part.charge),
                         Basic(lowercase(replace(part.mass.name,"ZERO"=>"0"))), 
                         Basic(lowercase(replace(part.width.name,"ZERO"=>"0"))) )
    push!( particle_list, new_part )
    push!( particle_name_dict, to_qgraf_name(part.name) => new_part )
    push!( particle_kf_dict, part.pdg_code => new_part )
  end # for part


  interaction_list = Vector{Interaction}()
  sorted_kf_list_dict = Dict{Vector{Int64},Interaction}()
  for vert in py_model.all_vertices
    # UFO convention is outgoing particle, but in our convention is incoming
    vert_link_list = map( p_ -> particle_name_dict[to_qgraf_name(p_.antiname)], vert.particles )
    
    has_yt = check_has_yt( vert.couplings )
    QCD_order, QED_order, SPC_order = extract_QCD_QED_order( vert.couplings )
    CTcoeff = calculate_CTcoeff( vert_link_list, has_yt, QCD_order )

    color_row_list = map( Basic, vert.color )
    color_dim = length(color_row_list)

    lorentz_col_list = map( l_ -> Basic(l_.structure), vert.lorentz )
    lorentz_dim = length(lorentz_col_list)

    couplings_matrix = extract_couplings_matrix( vert.couplings, color_dim, lorentz_dim, CTcoeff )

    new_interaction = Interaction( vert.name, vert_link_list, 
                                   color_row_list, lorentz_col_list, couplings_matrix,
                                   QCD_order, QED_order, SPC_order )

    push!( interaction_list, new_interaction )

    sorted_kf_list = sort( map( p_ -> p_.kf, vert_link_list ) )
    push!( sorted_kf_list_dict, sorted_kf_list => new_interaction )
  end # for vert

  # establish line counter-terms
  for part in particle_list
    if is_not_colorful(part) || part.kf < 0 
      continue
    end # if

    vert_name = part.name*"CT"
    vert_link_list = [ part, particle_name_dict[part.antiname] ]

    color_row_list, lorentz_col_list, couplings_matrix = generate_color_lorentz_couplings( part )

    new_interaction = Interaction( vert_name, vert_link_list, 
                                   color_row_list, lorentz_col_list, couplings_matrix,
                                   0 #=QCD_order=#, 0#=QED_order=#, 0#=SPC_order=# )

    push!( interaction_list, new_interaction )

    sorted_kf_list = sort( map( p_ -> p_.kf, vert_link_list ) )
    push!( sorted_kf_list_dict, sorted_kf_list => new_interaction )
  end # for part


  parameter_dict = Dict{Basic,Basic}() 
  for param in py_model.all_parameters
    param_name = lowercase( param.name )
    if param_name == "zero" 
      continue
    end # if
    param_value = lowercase(string( param.value ))

    param_value = replace( param_value, "complex(0,1)" => "I" )
    param_value = replace( param_value, "cmath.cos" => "cos" )
    param_value = replace( param_value, "cmath.sin" => "sin" )
    param_value = replace( param_value, "cmath.sqrt(2)" => "sqrt2" )
    param_value = replace( param_value, "cmath.sqrt" => "sqrt" )
    param_value = replace( param_value, "cmath.pi" => "pi" )
    param_value = replace( param_value, "complexconjugate" => "conj" )
    param_value = replace( param_value, "**" => "^" )

    push!( parameter_dict, Basic(param_name) => Basic(param_value) )
  end # for param

  coupling_dict = Dict{Basic, Basic}()
  for coupling ∈ py_model.all_couplings
    coupling_name = lowercase( coupling.name )
    coupling_value = lowercase( coupling.value )

    coupling_value = replace( coupling_value, "complex(0,1)" => "I" )
    coupling_value = replace( coupling_value, "cmath.cos" => "cos" )
    coupling_value = replace( coupling_value, "cmath.sin" => "sin" )
    coupling_value = replace( coupling_value, "cmath.sqrt(2)" => "sqrt2" )
    coupling_value = replace( coupling_value, "cmath.sqrt" => "sqrt" )
    coupling_value = replace( coupling_value, "cmath.pi" => "pi" )
    coupling_value = replace( coupling_value, "complexconjugate" => "conj" )
    coupling_value = replace( coupling_value, "**" => "^" )

    coupling_dict[Basic(coupling_name)] = Basic(coupling_value)
  end # for coupling

  #--------------------------------------
  # Universe Model instance
  model = Model( model_name, input["unitary_gauge"]::Bool, 
                 particle_list, particle_name_dict, particle_kf_dict, 
                 interaction_list, sorted_kf_list_dict, parameter_dict, coupling_dict )
  #--------------------------------------

  return model

end # function readin_model





