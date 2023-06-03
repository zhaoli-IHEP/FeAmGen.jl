
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


