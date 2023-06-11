
using FeAmGen

###########################
function main()::Nothing
###########################

  dir_name = "t_TO_t_4Loop/t_TO_t_4Loop_visuals"
  tex_list = filter( x->endswith(x,"tex"), readdir(dir_name) )
  println( "Show the non-planar diagrams" )
  for one_tex in tex_list
    if !is_planar( joinpath(dir_name,one_tex) )
      println( one_tex )
    end # if
  end # for one_tex

  return nothing

end # function main

########
main()
########
