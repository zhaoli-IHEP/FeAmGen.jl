using FeynUtils
using SymEngine

#########################
function main()::Nothing
#########################

  dir_name = "Wplus_t_TO_Wplus_t_3Loop/Wplus_t_TO_Wplus_t_3Loop_visuals"
  tex_list = filter( s->endswith(s,"tex"), readdir(dir_name) )

  the_index_list = Vector{Basic}()
  for one_tex in tex_list
    file = open( joinpath(dir_name,one_tex), "r" )
    content = read( file, String )
    close( file )

    line_list = filter( x->!isempty(x)&&x[1]=='v', string.(split( content, "\n" )) )
    bottom_list = filter( x->findfirst("label' = \\(b\\)",x)!=nothing, line_list )
    internal_bottom_list = filter( x->findfirst("q",x)==nothing, bottom_list )
    if !isempty(internal_bottom_list)
      index = Basic( one_tex[15:end][1:(end-4)] ) 
      println( "Has internal_bottom for cut: #$(index)" )
      push!( the_index_list, index )
    end # if
  end # for index
  sort!( the_index_list )

  @show the_index_list length(the_index_list)
  println()

  up_index_list = Vector{Basic}()
  for index in the_index_list
    file = open( joinpath( dir_name, "visual_diagram$(index).tex" ), "r" )
    content = read( file, String )
    close( file )

    line_list = filter( x->!isempty(x)&&x[1]=='v', string.(split( content, "\n" )) )
    up_list = filter( x->findfirst("label' = \\(u\\)",x)!=nothing, line_list )

    if !isempty(up_list)
      println( "Has up-quark: #$(index)" )
      push!( up_index_list, index )
    end # if

  end # for index

  @show up_index_list length(up_index_list)
  println()

  return nothing

end # function main


########
main()
########
