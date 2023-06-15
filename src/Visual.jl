
#####################################################
"""
    get_line_style_str( part::Particle )::String `(used only in generate_visual_graph)`

Explain the particle as a string in the relevant propagator.
"""
function get_line_style_str( part::Particle )::String
#####################################################

  if part.spin == :scalar 
    if is_not_majorana(part) && part.kf > 0 
      return "charged scalar"
    elseif is_not_majorana(part) && part.kf < 0 
      return "anti charged scalar"
    else
      return "scalar"
    end # if
  elseif part.spin == :fermion 
    if is_not_majorana(part) && part.kf > 0 
      return "fermion"
    elseif is_not_majorana(part) && part.kf < 0 
      return "anti fermion"
    else
      return "majorana"
    end # if
  elseif part.spin == :vector
    if is_gluon(part) 
      return "gluon"
    else
      return "boson"
    end # if
  elseif part.spin == :ghost 
    return "ghost"
  else
    error( "Spin exception!" )
  end # if

end # function get_line_style_str



#######################################################################
"""
    generate_visual_graph( g::Graph, model::Model )::String

Generate the string representing the relevant graph `g` in the tizk-feynamn syntax.
"""
function generate_visual_graph( g::Graph, model::Model )::String
#######################################################################

  result_str = 
    "\\begin{figure}[!htb] \n"*
    "\\begin{center} \n"*
    "\\feynmandiagram [large] { \n"*
    "  %$(g) \n"

  for att in g.property
    result_str *= 
    "    %$(att[1]) => $(att[2]) \n"
  end # for att

  for edge in g.edge_list
    result_str *= 
    "  %$(string(edge)) \n"
    for att in edge.property
      result_str *= 
      "    %$(att[1]) => $(att[2]) \n"
    end # for att
  end # for edge

  for vert in g.node_list
    result_str *= 
    "  %$(string(vert)) \n"
    for att in vert.property
      result_str *= 
      "    %$(att[1]) => $(att[2]) \n"
    end # for att
  end # for vert

  QCDct_str_dict = Dict( 0 => "", 1 => " [crossed dot]", 2 => " [red, label = DOUBLE, crossed dot]" )

  edge_list = g.edge_list
  for edge in edge_list
    src_mark = edge.src_node_mark
    dst_mark = edge.dst_node_mark
    mark_pair_set = Set([src_mark,dst_mark])

    parallel_edge_list = filter( e_ -> Set([e_.src_node_mark,e_.dst_node_mark]) == mark_pair_set, edge_list )

    if length(parallel_edge_list) <= 1
      half_circle_option = ""
    else 
      parallel_edge_mark_list = map( e_ -> e_.property[:mark], parallel_edge_list )   
      max_mark = max(parallel_edge_mark_list...)
      min_mark = min(parallel_edge_mark_list...)
      if edge.property[:mark] == max_mark
        half_circle_option = src_mark > dst_mark ? ", half left" : ", half right"
      elseif edge.property[:mark] == min_mark
        half_circle_option = src_mark > dst_mark ? ", half right" : ", half left"
      else # in the middle
        half_circle_option = "" 
      end # if
    end # end if

    src_QCDct_str = QCDct_str_dict[ get_node_mark_prop( g, src_mark, :QCDct_order ) ]

    dst_QCDct_str = QCDct_str_dict[ get_node_mark_prop( g, dst_mark, :QCDct_order ) ]

    edge_style_str = get_line_style_str(edge.property[:particle])

    tadpole_option = ""
    if src_mark == dst_mark 
      tadpole_option = ", loop, min distance=3cm"
    end # if
   
    src_external_marking_str = ""
    if get_node_mark_prop( g, src_mark, :style ) == "External"
      src_external_marking_str = " [particle = \\($(src_mark)\\)]"
    end # if
    dst_external_marking_str = ""
    if get_node_mark_prop( g, dst_mark, :style ) == "External"
      dst_external_marking_str = " [particle = \\($(dst_mark)\\)]"
    end # if


    particle_name = edge.property[:particle].name
    particle_name = replace( particle_name, "plus" => "^{+}" )
    particle_name = replace( particle_name, "minus" => "^{-}" )
    particle_name = replace( particle_name, "ve" => "\\nu_{e}" )
    particle_name = replace( particle_name, "mu" => "\\mu" )
    particle_name = replace( particle_name, "ta" => "\\tau" )
    particle_name = replace( particle_name, "vm" => "\\nu_{\\mu}" )
    particle_name = replace( particle_name, "vt" => "\\nu_{\\tau}" )
    particle_name = replace( particle_name, "^a" => "\\gamma" )
    if length(particle_name) > 3 && particle_name[end-2:end] == "bar"
      particle_name = "\\overline{"*particle_name[1:end-3]*"}"
    end # if

    mom_str = replace( string(edge.property[:momentum]), r"([Kkq]+)(\d+)" => s"\1_{\2}" )

    result_str *= 
      "$(get_node_mark_prop(g,src_mark,:name))$(src_QCDct_str)$(src_external_marking_str) -- [$(edge_style_str)$(half_circle_option)$(tadpole_option), edge label' = \\($(particle_name)\\), momentum = \\($(mom_str)\\) ] $(get_node_mark_prop(g,dst_mark,:name))$(dst_QCDct_str)$(dst_external_marking_str), \n"
  end # for edge

  result_str *= """
  };
  \\end{center}
  \\caption{Diagram$(g.property[:diagram_index]), Sign: $(g.property[:sign]), Symmetry factor: $(g.property[:symmetry_factor])}
  \\end{figure}
  \\newpage

  """

  ## remove redundant label for `[crossed dot]`
  v_crossed_dot_regex = r"v[1-9]\d* \[crossed dot\]"
  findnext_start_pointer = 1
  v_crossed_dot_position = findnext( v_crossed_dot_regex, result_str, findnext_start_pointer )
  crossed_dot_index_list = Int[]
  while !isnothing( v_crossed_dot_position )
    index = parse( Int, match( r"[1-9]\d*", result_str[v_crossed_dot_position] ).match )
    if index ∈ crossed_dot_index_list
      result_str = result_str[begin:(first(v_crossed_dot_position)-1)] *
                    replace( result_str[v_crossed_dot_position], " [crossed dot]" => "") *
                    result_str[(last(v_crossed_dot_position)+1):end]
    else
      push!( crossed_dot_index_list, index )
    end # if
    findnext_start_pointer = first( v_crossed_dot_position ) + 1
    v_crossed_dot_position = findnext( v_crossed_dot_regex, result_str, findnext_start_pointer )
  end # while
  @assert all( index->(length∘findall)("v$index [crossed dot]",result_str)==1, crossed_dot_index_list )
  # num_crossed_dot_str = (length∘findall)( "[crossed dot]", result_str )
  # first_range = findfirst( "[crossed dot]", result_str )
  # while num_crossed_dot_str > 1
  #   second_range = findnext( "[crossed dot]", result_str, last(first_range) )
  #   result_str = result_str[begin:first(second_range)-1] * result_str[last(second_range)+1:end]
  #   num_crossed_dot_str -= 1
  # end # while

  return result_str

end # function generate_visual_graph

