
#################################################
"""
    electron_num( part::Particle )::Int64

Count the electron number (lepton number in 1st generation ).
"""
function electron_num( part::Particle )::Int64
  if part.kf in [11,12]
    return 1
  end
  if part.kf in [-11,-12]
    return -1
  end
  return 0
end # function electron_num
#################################################
"""
    mu_num( part::Particle )::Int64

Count the muon number (lepton number in 2nd generation ).
"""
function mu_num( part::Particle )::Int64
  if part.kf in [13,14]
    return 1
  end
  if part.kf in [-13,-14]
    return -1
  end
  return 0
end # function mu_num
#################################################
"""
    tau_num( part::Particle )::Int64

Count the tauon number (lepton number in 3rd generation ).
"""
function tau_num( part::Particle )::Int64
  if part.kf in [15,16]
    return 1
  end
  if part.kf in [-15,-16]
    return -1
  end
  return 0
end # function tau_num
#################################################


#################################################
"""
    quark1st_num( part::Particle )::Int64

Count the quark flavor number in 1st generation.
This is used for the case that has no flavor violation.
"""
function quark1st_num( part::Particle )::Int64
  if part.kf in [1,2]
    return 1
  end
  if part.kf in [-1,-2]
    return -1
  end
  return 0
end # function quark1st_num
#################################################
"""
    quark2nd_num( part::Particle )::Int64

Count the quark flavor number in 2nd generation.
This is used for the case that has no flavor violation.
"""
function quark2nd_num( part::Particle )::Int64
  if part.kf in [3,4]
    return 1
  end
  if part.kf in [-3,-4]
    return -1
  end
  return 0
end # function quark2nd_num
#################################################
"""
    quark3rd_num( part::Particle )::Int64

Count the quark flavor number in 3rd generation.
This is used for the case that has no flavor violation.
"""
function quark3rd_num( part::Particle )::Int64
  if part.kf in [5,6]
    return 1
  end
  if part.kf in [-5,-6]
    return -1
  end
  return 0
end # function quark3rd_num
#################################################



#################################################
"""
    filter_lepton_generations( proc_str::Union{String,Nothing}, n_inc::Int64, part_dict::Dict{String,Particle} )::Union{String,Nothing}

This function is used when the option `AllowLeptonNumberViolation` is turned on.
It will check if the process represented by the input `proc_str` has violation the lepton number.
"""
function filter_lepton_generations( 
    proc_str::Union{String,Nothing}, n_inc::Int64, 
    part_dict::Dict{String,Particle} )::Union{String,Nothing}
#################################################
  if proc_str == nothing
    return nothing
  end

  part_str_list = split( proc_str, "," )
  part_list = map( s_ -> part_dict[s_], part_str_list )

  inc_electron_num_list = map( p_ -> electron_num(p_), part_list[1:n_inc] )
  inc_electron_num_sum = sum( inc_electron_num_list ) 
  out_electron_num_list = map( p_ -> electron_num(p_), part_list[n_inc+1:end] )
  out_electron_num_sum = sum( out_electron_num_list ) 
  if inc_electron_num_sum != out_electron_num_sum
    return nothing
  end

  inc_mu_num_list = map( p_ -> mu_num(p_), part_list[1:n_inc] )
  inc_mu_num_sum = sum( inc_mu_num_list ) 
  out_mu_num_list = map( p_ -> mu_num(p_), part_list[n_inc+1:end] )
  out_mu_num_sum = sum( out_mu_num_list ) 
  if inc_mu_num_sum != out_mu_num_sum
    return nothing
  end

  inc_tau_num_list = map( p_ -> tau_num(p_), part_list[1:n_inc] )
  inc_tau_num_sum = sum( inc_tau_num_list ) 
  out_tau_num_list = map( p_ -> tau_num(p_), part_list[n_inc+1:end] )
  out_tau_num_sum = sum( out_tau_num_list ) 
  if inc_tau_num_sum != out_tau_num_sum
    return nothing
  end

  return proc_str
end #function filter_lepton_generations






#################################################
"""
    filter_quark_generations( proc_str::Union{String,Nothing}, n_inc::Int64, part_dict::Dict{String,Particle} )::Union{String,Nothing}

This function is used when the option `AllowQuarkGenerationViolation` is turned on.
It will check if the process represented by the input `proc_str` has violation the quark flavor.
"""
function filter_quark_generations( 
    proc_str::Union{String,Nothing}, n_inc::Int64, 
    part_dict::Dict{String,Particle} )::Union{String,Nothing}
#################################################
  if proc_str == nothing
    return nothing
  end

  part_str_list = split( proc_str, "," )
  part_list = map( s_ -> part_dict[s_], part_str_list )

  inc_quark1st_num_list = map( p_ -> quark1st_num(p_), part_list[1:n_inc] )
  inc_quark1st_num_sum = sum( inc_quark1st_num_list ) 
  out_quark1st_num_list = map( p_ -> quark1st_num(p_), part_list[n_inc+1:end] )
  out_quark1st_num_sum = sum( out_quark1st_num_list ) 
  if inc_quark1st_num_sum != out_quark1st_num_sum
    return nothing
  end

  inc_quark2nd_num_list = map( p_ -> quark2nd_num(p_), part_list[1:n_inc] )
  inc_quark2nd_num_sum = sum( inc_quark2nd_num_list ) 
  out_quark2nd_num_list = map( p_ -> quark2nd_num(p_), part_list[n_inc+1:end] )
  out_quark2nd_num_sum = sum( out_quark2nd_num_list ) 
  if inc_quark2nd_num_sum != out_quark2nd_num_sum
    return nothing
  end

  inc_quark3rd_num_list = map( p_ -> quark3rd_num(p_), part_list[1:n_inc] )
  inc_quark3rd_num_sum = sum( inc_quark3rd_num_list ) 
  out_quark3rd_num_list = map( p_ -> quark3rd_num(p_), part_list[n_inc+1:end] )
  out_quark3rd_num_sum = sum( out_quark3rd_num_list ) 
  if inc_quark3rd_num_sum != out_quark3rd_num_sum
    return nothing
  end

  return proc_str
end #function filter_quark_generations



#################################################
"""
    filter_charge( proc_str::Union{String,Nothing}, n_inc::Int64, part_dict::Dict{String,Particle} )::Union{String,Nothing}

This function is used to check if the process represented by the input `proc_str` has violation the charge conservation.
"""
function filter_charge( 
    proc_str::Union{String,Nothing}, n_inc::Int64, 
    part_dict::Dict{String,Particle} )::Union{String,Nothing}
#################################################
  if proc_str == nothing
    return nothing
  end

  part_str_list = split( proc_str, "," )
  part_list = map( s_ -> part_dict[s_], part_str_list )

  inc_charge_list = map( p_ -> p_.charge, part_list[1:n_inc] )
  inc_charge_sum = sum( inc_charge_list ) 
  out_charge_list = map( p_ -> p_.charge, part_list[n_inc+1:end] )
  out_charge_sum = sum( out_charge_list ) 
  if inc_charge_sum != out_charge_sum
    return nothing
  end

  return proc_str
end #function filter_charge









#############################################################
"""
    simple_readin_model( model_name::String )::Dict{String,Particle}

Read-in the model detail from python model file. For seed program, we only need particle list.
The directory of model files are supposed in ".".
"""
function simple_readin_model(
  model_name::String;
  model_paths::Vector{String}=[pwd()]
)::Dict{String,Particle}
############################################################################################

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
  model = pyimport( "$(model_name).object_library" )

  particle_dict = Dict{String,Particle}()
  for py_part in model.all_particles
    new_part = Particle( py_part.pdg_code, to_qgraf_name(py_part.name), to_qgraf_name(py_part.antiname), 
                         spin_dict[py_part.spin], color_dict[py_part.color], charge_convert(py_part.charge), 
                         Basic(replace(py_part.mass.name,"ZERO"=>"0")),
                         Basic(replace(py_part.width.name,"ZERO"=>"0")) )
    push!( particle_dict, to_qgraf_name(py_part.name) => new_part )
  end # for part

  return particle_dict

end # function simple_readin_model





# ###########################################
# """
#     tensor_product( aa::Vector{String}, bb::Vector{String}, xx::Vector{String}... )::Vector{String} 

# Generalize tensor_product to multiple string lists/arrays.
# """
# function tensor_product( aa::Vector{String}, bb::Vector{String}, xx::Vector{String}... )::Vector{String} 
# ###########################################
#   return tensor_product( tensor_product( aa, bb ), xx... )
# end # function tensor_product



# ###########################################
# """
#     tensor_product( str_list1::Vector{String}, str_list2::Vector{String} )::Vector{String}

# Inside tensor_product( x, y ), x and y are two string lists/arrays.
# This function calculate the tensor production of two arrays.

# # Examples
# ```julia
# @testset "tensor_product" begin
#   target_in = Array{String}[ String["a","b","c"], String["d","e","f"] ]
#   target_out = String["a,d","a,e","a,f","b,d","b,e","b,f","c,d","c,e","c,f"]
#   @test FeAmGen.tensor_product( target_in... ) == target_out
#   target_in = Array{String}[ String["a","b"], String["c","d"], String["e","f"] ]
#   target_out = String["a,c,e","a,c,f","a,d,e","a,d,f","b,c,e","b,c,f","b,d,e","b,d,f"]
#   @test FeAmGen.tensor_product( target_in... ) == target_out
# end # @testset
# ```
# """
# function tensor_product( 
#     str_list1::Vector{String}, 
#     str_list2::Vector{String}
# )::Vector{String}
# ###########################################
#   res = Vector{String}()
#   for str1 in str_list1, str2 in str_list2
#     push!(res,str1*","*str2)
#   end
#   return res
# end # function tensor_product







###########################################
"""
    expand_parton( inc_str_list::Vector{String}, out_str_list::Vector{String}, parton_str_list::Vector{String})::Vector{String}

Sometimes we may use `parton` as one of the particles to represent the processes. 
This function can expand `parton` to generate several specific processes.
However, this expansion results are not checked for availability.

# Examples
```julia
@testset "expand_parton" begin
  @test FeAmGen.expand_parton( String["u","parton"], String["d","g"], String["u","d","g"] ) == String["u,u,d,g", "u,d,d,g", "u,g,d,g"]
end # @testset
```
"""
function expand_parton( 
    inc_str_list::Vector{String}, 
    out_str_list::Vector{String},
    parton_str_list::Vector{String}
)::Vector{String}
###########################################

  # suppose the list is [ "u", "parton", "d", "g" ]
  incout_str_list = [ inc_str_list; out_str_list ]

  # now explain parton into [ ["u"], ["u","d","g"], ["d"], ["g"] ]
  incout_parton_str_list = map( s_ -> s_ == "parton" ? parton_str_list : [s_], incout_str_list )

  # generate the tensor product
  # product_list = tensor_product( incout_parton_str_list... )
  product_list = vectorized_tensor_product_String( incout_parton_str_list... )

  return product_list
end # function expand_parton





################################################################
"""
    sort_proc_str( proc_str::Union{String,Nothing}, n_inc::Int64 )::Union{String,Nothing}

Sort the process string in the lexicographic ordering for the incoming part and outgoing part respectively.

# Examples
```julia
@testset "sort_proc_str" begin
  @test FeAmGen.sort_proc_str( "u,d,g,d,u", 2 ) == "d,u,d,g,u"
end # @testset
```
"""
function sort_proc_str( proc_str::Union{String,Nothing}, n_inc::Int64 )::Union{String,Nothing}
################################################################
  if proc_str == nothing
    return nothing
  end
  part_str_list = split( proc_str, "," )
  inc_str_list = convert( Vector{String}, part_str_list[1:n_inc] )
  out_str_list = convert( Vector{String}, part_str_list[n_inc+1:end] )
  sort!( inc_str_list )
  sort!( out_str_list )

  return join( [ inc_str_list; out_str_list ], "," )
end # function sort_proc_str




########################################################################
"""
    get_list_quoted_str( str_list::Vector{String} )::String

Convert from a list of string into a string of quoted strings.

# Examples 
```julia
@testset "get_list_quoted_str" begin
  @test FeAmGen.get_list_quoted_str( String["a","b","c"] ) == "[ \"a\",\"b\",\"c\" ]"
end # @testset
```
"""
function get_list_quoted_str( str_list::Vector{String} )::String
########################################################################

  quoted_str_list = map( s_ -> "\""*s_*"\"", str_list )

  return "[ "*join( quoted_str_list, ",")*" ]"
end # function get_list_quoted_str


####################################################################################
"""
    write_card( parton_proc_str::String, proc_str::String, n_inc::Int64, input::Dict{Any,Any} )::String

Generate the YAML input card for the specific process.
"""
function write_card( parton_proc_str::String, proc_str::String, n_inc::Int64, input::Dict{Any,Any} )::String
####################################################################################

  part_str_list = split( proc_str, "," )
  proc_name = join( [ part_str_list[1:n_inc]; "TO"; part_str_list[n_inc+1:end] ], "_" )

  converted_incomings = convert( Vector{String}, part_str_list[1:n_inc] )
  converted_outgoings = convert( Vector{String}, part_str_list[n_inc+1:end] )
  inc_part_list_str = get_list_quoted_str( converted_incomings )
  out_part_list_str = get_list_quoted_str( converted_outgoings )

  proc_card_path = joinpath( parton_proc_str, "$(proc_name).yaml" )
  file = open( proc_card_path, "w+" )
  write( file, """
  # input file generated by seed program
  
  model_name: \"$(input["model_name"])\"
  
  # use unitary_gauge for internal vector boson
  unitary_gauge: $(input["unitary_gauge"])
  
  # process information
  DropTadpole: $(input["DropTadpole"])              # drop tadpole?
  DropWFcorrection: $(input["DropWFcorrection"])         # drop WFcorrection?
  
  # number of loops
  n_loop: $(input["n_loop"])
  # order of QCD counter-term vertices
  QCDCT_order: $(input["QCDCT_order"])
  
  # order of QCD coupling gs in the amplitude
  Amp_QCD_order: $(input["Amp_QCD_order"])
  # order of QED coupling ee in the amplitude
  Amp_QED_order: $(input["Amp_QED_order"])
  # order of special coupling in the amplitude
  Amp_SPC_order: $(input["Amp_SPC_order"])
  
  # min ep power in the amplitude
  Amp_Min_Ep_Xpt: $(input["Amp_Min_Ep_Xpt"])
  # max ep power in the amplitude
  Amp_Max_Ep_Xpt: $(input["Amp_Max_Ep_Xpt"])
  
  # incoming and outgoing information
  incoming: $inc_part_list_str           # incoming particles
  outgoing: $out_part_list_str           # outgoing particles 
  
  # Used for simplifying amplitude expression, i.e. "result amplitude" = "true amplitude"/couplingfactor
  couplingfactor: \"1\"

  # Symmetry Configuration
  momentum_symmetry: $(isempty(input["momentum_symmetry"]) ? "[]" : input["momentum_symmetry"])
  color_symmetry: $(isempty(input["color_symmetry"]) ? "[]" : input["color_symmetry"])

  """ )
  close(file)

  return proc_card_path

end # function write_card










