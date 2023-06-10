
mutable struct Node
  mark::Int64
  property::Dict{Symbol,Any}
end # struct Node


mutable struct Edge
  mark::Int64
  src_node_mark::Int64
  dst_node_mark::Int64
  property::Dict{Symbol,Any}
end # struct Edge


mutable struct Graph
  node_list::Vector{Node}
  edge_list::Vector{Edge}
  property::Dict{Symbol,Any}
end # struct Graph


function graph()::Graph
  return Graph( Vector{Node}(), Vector{Edge}(), Dict{Symbol,Any}() )
end # function graph

function set_prop!( g::Graph, key::Symbol, val::Any )::Nothing
  g.property[key] = val 
  return nothing
end # function set_prop!

function set_prop!( g::Graph, new_property::Dict{Symbol,Any} )::Nothing
  g.property = (Dict∘union)( g.property, new_property ) 
  return nothing
end # function set_prop!

function get_prop( g::Graph, key::Symbol )::Any
  return g.property[key]  
end # function get_prop 



function add_node!( g::Graph, node::Node )::Nothing
  @assert (isempty∘filter)( n -> n.mark == node.mark, g.node_list )
  push!( g.node_list, node )
  return nothing
end # function add_node!

function add_node!( g::Graph, node_mark::Int64 )::Nothing
  @assert (isempty∘filter)( n -> n.mark == node_mark, g.node_list )
  push!( g.node_list, Node( node_mark, Dict{Symbol,Any}() ) )
  return nothing
end # function add_node!

function add_node!( g::Graph, node_mark::Int64, node_property::Dict{Symbol,Any} )::Nothing
  @assert (isempty∘filter)( n -> n.mark == node_mark, g.node_list )
  push!( g.node_list, Node( node_mark, node_property ) )
  return nothing
end # function add_node!



function add_edge!( g::Graph, edge::Edge )::Nothing 
  @assert (isempty∘filter)( e -> e.mark == edge.mark, g.edge_list )
  push!( g.edge_list, edge )
  return nothing
end # function add_edge!


function add_edge!( g::Graph, edge_mark::Int64, src_mark::Int64, dst_mark::Int64 )::Nothing
  @assert (isempty∘filter)( e -> e.mark == edge_mark, g.edge_list )
  push!( g.edge_list, Edge( edge_mark, src_mark, dst_mark, Dict{Symbol,Any}() ) )
  return nothing
end # function add_edge!

function add_edge!( g::Graph, edge_mark::Int64, src_mark::Int64, dst_mark::Int64, edge_property::Dict{Symbol,Any} )::Nothing
  @assert (isempty∘filter)( e -> e.mark == edge_mark, g.edge_list )
  push!( g.edge_list, Edge( edge_mark, src_mark, dst_mark, edge_property ) )
  return nothing
end # function add_edge!

function add_edge!( g::Graph, edge_mark::Int64, src_dst_mark_pair::Tuple{Int64,Int64}, edge_property::Dict{Symbol,Any} )::Nothing
  add_edge!( g, edge_mark, src_dst_mark_pair[1], src_dst_mark_pair[2], edge_property )
  return nothing
end # function add_edge!





function n_node( g::Graph )::Int64
  return length(g.node_list)
end # function n_node

function n_edge( g::Graph )::Int64
  return length(g.edge_list)
end # function n_edge



function get_node_mark_prop( g::Graph, node_mark::Int64, key::Symbol )::Any
  for node in g.node_list
    if node.mark == node_mark
      return node.property[key]
    end # if
  end # for node
  error("Cannot find property $(key)")
  return 0
end # function get_node_mark_prop


function get_node_index_prop( g::Graph, node_index::Int64, key::Symbol )::Any
  return g.node_list[node_index].property[key]
end # function get_node_index_prop


function get_edge_mark_prop( g::Graph, edge_mark::Int64, key::Symbol )::Any
  for edge in g.edge_list
    if edge.mark == edge_mark
      return edge.property[key]
    end # if
  end # for node
  error("Cannot find property $(key)")
  return 0
end # function get_edge_mark_prop


function get_edge_index_prop( g::Graph, edge_index::Int64, key::Symbol )::Any
  return g.edge_list[edge_index].property[key]
end # function get_edge_index_prop




function get_in_edge_list( g::Graph, node_index::Int64 )::Vector{Edge}

  node_mark = g.node_list[node_index].mark

  return filter( e_ -> e_.dst_node_mark == node_mark, g.edge_list )

end # function get_in_edge_list

function get_out_edge_list( g::Graph, node_index::Int64 )::Vector{Edge}

  node_mark = g.node_list[node_index].mark

  return filter( e_ -> e_.src_node_mark == node_mark, g.edge_list )

end # function get_out_edge_list


###############################################################
@inline is_external( edge::Edge)::Bool = edge.property[:style] == "External"
###############################################################


###############################################################
# Created by Quan-feng WU
# 2023-06-10
function is_planar( visual_tex_file::String )::Bool
###############################################################
  @assert isfile( visual_tex_file )
  contents = readlines( visual_tex_file )
  filter!( startswith(r"v[1-9]\d*"), contents )

  edges_list = Tuple{Int,Int}[]
  for one_line ∈ contents
    v_range_list = findall( r"v[1-9]\d*", one_line )
    vertex_index_list = [ parse( Int, one_line[v_range][2:end] )
                            for v_range ∈ v_range_list ]
    n_vertices = length( vertex_index_list )
    for ii ∈ 1:n_vertices-1
      vi = vertex_index_list[ii]
      vj = vertex_index_list[ii+1]
      if vi ≤ vj
        push!( edges_list, (vi, vj) )
      else
        push!( edges_list, (vj, vi) )
      end # if
    end # for ii
  end # for one_line

  vertices_list = (sort∘union)( edges_list... )
  vertices_count = [ sum( e->count(==(vi),e), edges_list )
                      for vi ∈ vertices_list ]
  v_inf = maximum( vertices_list ) + 1
  v_ext_indices = findall( ==(1), vertices_count )
  for v_ext ∈ vertices_list[v_ext_indices]
    push!( edges_list, (v_ext, v_inf) )
  end # for v_ext
  push!( vertices_list, v_inf )

  n_vertices = length( vertices_list )
  for (ee, one_edge) ∈ enumerate( edges_list )
    vi, vj = one_edge
    vi = findfirst( ==(vi), vertices_list )
    vj = findfirst( ==(vj), vertices_list )
    edges_list[ee] = (vi, vj)
  end # for (ee, one_edge)

  next_vertex = n_vertices + 1
  to_be_added_edges_list = Tuple{Int,Int}[]
  self_loop_edge_indices = findall( iijj->first(iijj)==last(iijj), edges_list )
  for ee ∈ self_loop_edge_indices
    vi = first( edges_list[ee] )
    push!( to_be_added_edges_list,
            (vi, next_vertex),
            (next_vertex, next_vertex+1),
            (vi, next_vertex+1) )
    next_vertex += 2
    n_vertices += 2
  end # for ee
  deleteat!( edges_list, self_loop_edge_indices )
  edges_list = vcat( edges_list, to_be_added_edges_list )

  unique_edges_list = unique( edges_list )
  to_be_added_edges_list = Tuple{Int,Int}[]
  for unique_edge ∈ unique_edges_list
    vi, vj = unique_edge
    duplicate_number = count( ==(unique_edge), edges_list )
    @assert duplicate_number ≥ 1
    while duplicate_number > 1
      push!( to_be_added_edges_list, 
              (vi, next_vertex),
              (vj, next_vertex) )
      next_vertex += 1
      n_vertices += 1
      duplicate_number -= 1
    end # while
  end # for unique_edge
  edges_list = vcat( unique_edges_list, to_be_added_edges_list )

  adj_mat = zeros( Int, n_vertices, n_vertices )
  for one_edge ∈ edges_list
    ii, jj = one_edge
    adj_mat[ii, jj] += 1
    adj_mat[jj, ii] += 1
  end # for one_edge
  @assert all( ele->iszero(ele)||ele==1, adj_mat )
  @assert all( iszero, adj_mat[ii, ii] for ii ∈ 1:n_vertices )

  adj_mat_str = "n=$n_vertices\n" * join( [ join(row) for row ∈ eachrow(adj_mat) ], "\n" )

  write( "adj_mat.txt", adj_mat_str )

  amtog_cmd = pipeline( `$(amtog()) adj_mat.txt input.g6`, stdout="amtog.log" )
  run( amtog_cmd )

  planarg_cmd = pipeline( `$(planarg()) input.g6 planar_result.txt` )
  run( planarg_cmd )

  planar_list = readlines( "planar_result.txt" )
  isempty( planar_list ) && return false
  return true

end # function is_planar
