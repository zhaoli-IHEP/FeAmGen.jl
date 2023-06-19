
using FeAmGen

###########################
function main()::Nothing
###########################

  dir_name = "t_TO_t_4Loop/t_TO_t_4Loop_visuals"
  tex_list = filter( x->endswith(x,"tex"), readdir(dir_name) )
  println( "Show the non-planar diagrams" )
  nonplanar_list = Vector{String}()
  for one_tex in tex_list
    if !is_planar( joinpath(dir_name,one_tex) )
      push!( nonplanar_list, one_tex )
    end # if
  end # for one_tex

  nonplanar_index_list = map( (string∘first), [[ m.match for m in eachmatch(r"\d+", str)] for str in nonplanar_list ] )
  @show nonplanar_index_list 


  return nothing

end # function main

########
main()
########
