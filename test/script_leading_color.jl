using JLD2
using SymEngine
using FeynUtils

##########################
function main()::Nothing
##########################

  @vars nc ca cf 
  dict = Dict{Basic,Basic}([ ca => nc, cf => (nc^2-1)/(2*nc) ])

  dir_name = "t_TO_t_4Loop/t_TO_t_4Loop_amplitudes"
  jld_list = filter( x->endswith(x,"jld2"), readdir(dir_name) )
  println( "Show the diagrams containing highest power of Nc, i.e. leading color." )
  max_xpt = zero(Int64)
  for one_jld in jld_list
    file = jldopen( joinpath(dir_name,one_jld), "r" )
    color_list = (to_Basic∘read)( file, "amp_color_list" )
    close( file )
    color_list = map( x->(expand∘subs)( x, dict ), color_list )
    for one_color in color_list
      the_max_xpt = get_exponent(one_color,nc,:Max)
      max_xpt = max( max_xpt, the_max_xpt )
    end # for one_color
  end # for one_jld

  println( "Leading color is Nc^$(max_xpt)" )

  leading_color_list = Vector{String}()
  for one_jld in jld_list
    file = jldopen( joinpath(dir_name,one_jld), "r" )
    color_list = (to_Basic∘read)( file, "amp_color_list" )
    close( file )
    color_list = map( x->(expand∘subs)( x, dict ), color_list )
    for one_color in color_list
      the_coeff = coeff( one_color, nc, Basic(max_xpt) ) 
      if !iszero(the_coeff)
        push!( leading_color_list, (first∘splitext)(one_jld) )
        break
      end # if
    end # for one_color
  end # for one_jld

  leading_color_index_list = map( (string∘first), [[ m.match for m in eachmatch(r"\d+", str)] for str in leading_color_list ] )
  @show leading_color_index_list 

  return nothing

end # function main


##############
main()
##############
