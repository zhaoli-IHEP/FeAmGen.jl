using FeynUtils
using SymEngine

#########################
function main()::Nothing
#########################

  root, dirs, files = (first∘collect∘walkdir)("./Wplus_t_TO_Wplus_t_3Loop_visuals")
  tex_list = filter( s->endswith(s,".tex"), files )

  the_index_list = Vector{Basic}()
  for one_tex in tex_list
    file_name = "Wplus_t_TO_Wplus_t_3Loop_visuals/$(one_tex)"
    file = open( file_name, "r" )
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


  up_index_list = Vector{Basic}()
  for index in the_index_list
    file_name = "Wplus_t_TO_Wplus_t_3Loop_visuals/visual_diagram$(index).tex"
    file = open( file_name, "r" )
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

  bk_mkdir( "Wplus_t_TO_Wplus_t_3Loop_amplitudes_nf" )
  for index in up_index_list
    cp( "Wplus_t_TO_Wplus_t_3Loop_amplitudes/amplitude_diagram$(index).out", 
        "Wplus_t_TO_Wplus_t_3Loop_amplitudes_nf/amplitude_diagram$(index).out" ) 
    cp( "Wplus_t_TO_Wplus_t_3Loop_amplitudes/amplitude_diagram$(index).jld2",  
        "Wplus_t_TO_Wplus_t_3Loop_amplitudes_nf/amplitude_diagram$(index).jld2" ) 
  end # for index

  bk_mkdir( "Wplus_t_TO_Wplus_t_3Loop_visuals_nf" )
  cp( "Wplus_t_TO_Wplus_t_3Loop_visuals/generate_diagram_pdf.jl", 
      "Wplus_t_TO_Wplus_t_3Loop_visuals_nf/generate_diagram_pdf.jl" )
  cp( "Wplus_t_TO_Wplus_t_3Loop_visuals/tikz-feynman.sty", 
      "Wplus_t_TO_Wplus_t_3Loop_visuals_nf/tikz-feynman.sty" )
  for index in up_index_list
    cp( "Wplus_t_TO_Wplus_t_3Loop_visuals/visual_diagram$(index).tex", 
        "Wplus_t_TO_Wplus_t_3Loop_visuals_nf/visual_diagram$(index).tex" ) 
    cp( "Wplus_t_TO_Wplus_t_3Loop_visuals/expression_diagram$(index).out", 
        "Wplus_t_TO_Wplus_t_3Loop_visuals_nf/expression_diagram$(index).out" ) 
  end # for index


  return nothing

end # function main


########
main()
########
