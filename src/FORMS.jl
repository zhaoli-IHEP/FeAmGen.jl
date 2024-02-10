





# ##################################################################
# """
#     make_baseINC_script( graph::Graph )::String

# Prepare the baseINC script for containing kinematic relations.
# """
# function make_baseINC_script( graph::Graph )::String
# ###################################################################

#   result_str = string()

#   n_inc = graph.property[:n_inc]
#   n_out = graph.property[:n_out]
#   n_leg = n_inc+n_out

#   ext_edge_list = filter( is_external, graph.edge_list )
#   sorted_ext_edge_list = sort( ext_edge_list, by=x->x.property[:mark] ) 

#   #---------------------
#   if n_leg == 2

#     @assert sorted_ext_edge_list[1].property[:mark] == 1
#     mom1 = sorted_ext_edge_list[1].property[:momentum]
#     @assert sorted_ext_edge_list[2].property[:mark] == 2
#     mom2 = sorted_ext_edge_list[2].property[:momentum]
#     result_str = """
#     argument Levi;
#       id $(mom2) = $(mom1);
#     endargument;

#     """

#     return result_str
#   end # if
#   #---------------------


#   momN = sorted_ext_edge_list[n_leg].property[:momentum]
#   momNm1 = sorted_ext_edge_list[n_leg-1].property[:momentum]
#   momNm2 = sorted_ext_edge_list[n_leg-2].property[:momentum]

#   #-------------------------------------------------------------
#   sorted_notN_ext_edge_list = filter( x -> x.property[:mark] != n_leg, sorted_ext_edge_list )

#   result_str *= "id FV($(momN),rho?) = ";
#   for edge in sorted_notN_ext_edge_list
#     mom = edge.property[:momentum]
#     inc_sign = edge.property[:mark] <= n_inc ? 1 : (-1)
#     result_str *= "+($(inc_sign))*FV($(mom),rho)"
#   end # for edge
#   result_str *= ";\n"
 
#   result_str *= "id SP($(momN),rho?) = ";
#   for edge in sorted_notN_ext_edge_list
#     mom = edge.property[:momentum]
#     inc_sign = edge.property[:mark] <= n_inc ? 1 : (-1)
#     result_str *= "+($(inc_sign))*SP($(mom),rho)"
#   end # for edge
#   result_str *= ";\n"

#   result_str *= 
#     "id FermionChain( Spinor1?ILSPSET(int1?!{,$(n_leg)},mom1?,mass1?),"*
#     " ?vars1, GA($(momN)), ?vars2,"*
#     " Spinor2?IRSPSET(int2?!{,$(n_leg)},mom2?,mass2?) ) = \n"
#   for edge in sorted_notN_ext_edge_list
#     mom = edge.property[:momentum]
#     inc_sign = edge.property[:mark] <= n_inc ? 1 : (-1)
#     result_str *= 
#     "  +($(inc_sign))*FermionChain( Spinor1(int1,mom1,mass1), ?vars1, GA($(mom)), ?vars2, Spinor2(int2,mom2,mass2) )\n"
#   end # for edge
#   result_str *= ";\n"

#   #-------------------------------------------------------------
#   sorted_notNm1_ext_edge_list = filter( x -> x.property[:mark] != n_leg-1, sorted_ext_edge_list )
#   Nm1_sign = n_leg-1 <= n_inc ? (-1) : (+1)

#   result_str *=
#     "id FermionChain( Spinor?ILSPSET($(n_leg),mom?,mass?), ?vars1, GA($(momNm1)), ?vars2 ) = \n"
#   for edge in sorted_notNm1_ext_edge_list
#     mom = edge.property[:momentum]
#     inc_sign = edge.property[:mark] <= n_inc ? 1*Nm1_sign : (-1)*Nm1_sign
#     result_str *=
#     "  +($(inc_sign))*FermionChain( Spinor($(n_leg),mom,mass), ?vars1, GA($(mom)), ?vars2 )\n"
#   end # for edge
#   result_str *= ";\n"

#   result_str *=
#     "id FermionChain( ?vars1, GA($(momNm1)), ?vars2, Spinor?IRSPSET($(n_leg),mom?,mass?) ) = \n"
#   for edge in sorted_notNm1_ext_edge_list
#     mom = edge.property[:momentum]
#     inc_sign = edge.property[:mark] <= n_inc ? 1*Nm1_sign : (-1)*Nm1_sign
#     result_str *=
#     "  +($(inc_sign))*FermionChain( ?vars1, GA($(mom)), ?vars2, Spinor($(n_leg),mom,mass) )\n"
#   end # for edge
#   result_str *= ";\n"

#   #-------------------------------------------------------------
#   sorted_notNm2_ext_edge_list = filter( x -> x.property[:mark] != n_leg-2, sorted_ext_edge_list )
#   Nm2_sign = n_leg-2 <= n_inc ? (-1) : (+1)
  
#   result_str *=
#     "id FermionChain( Spinor1?ILSPSET($(n_leg-1),mom1?,mass1?), ?vars1, GA($(momNm2)), ?vars2, Spinor2?IRSPSET($(n_leg),mom2?,mass2?) ) = \n"
#   for edge in sorted_notNm2_ext_edge_list
#     mom = edge.property[:momentum]
#     inc_sign = edge.property[:mark] <= n_inc ? 1*Nm2_sign : (-1)*Nm2_sign
#     result_str *=
#     "  +($(inc_sign))*FermionChain( Spinor1($(n_leg-1),mom1,mass1), ?vars1, GA($(mom)), ?vars2, Spinor2($(n_leg),mom2,mass2) )\n"
#   end # for edge
#   result_str *= ";\n"

#   result_str *=
#     "id FermionChain( Spinor1?ILSPSET($(n_leg),mom1?,mass1?), ?vars1, GA($(momNm2)), ?vars2, Spinor2?IRSPSET($(n_leg-1),mom2?,mass2?) ) = \n"
#   for edge in sorted_notNm2_ext_edge_list
#     mom = edge.property[:momentum]
#     inc_sign = edge.property[:mark] <= n_inc ? 1*Nm2_sign : (-1)*Nm2_sign
#     result_str *=
#     "  +($(inc_sign))*FermionChain( Spinor1($(n_leg),mom1,mass1), ?vars1, GA($(mom)), ?vars2, Spinor2($(n_leg-1),mom2,mass2) )\n"
#   end # for edge
#   result_str *= ";\n"

#   #-----------------------------------------------------------------------------------
#   if sorted_ext_edge_list[n_leg].property[:particle].spin == :vector
#     result_str *=
#       "id FV($(momNm1),rho?)*VecEpsilon?{VecEp,VecEpC}($(n_leg),rho?,$(momN),r$(n_leg)?,mass?) = \n"
#     for index in 1:(n_leg-2)
#       edge = sorted_ext_edge_list[index]
#       mom = edge.property[:momentum]
#       inc_sign = index <= n_inc ? (+1) : (-1)
#       result_str *=
#       "  +($(inc_sign))*FV($(mom),rho)*VecEpsilon($(n_leg),rho,$(momN),r$(n_leg),mass)\n"
#     end # for index
#   end # if

#   result_str *= """
#     ; 
#     id FV(mom?,rho?)*VecEpsilon?{VecEp,VecEpC}(int?,rho?,mom?,mass?) = 0;
#     """

#   return result_str

# end # function make_baseINC_script









##############################################################################
# Created by Quan-feng WU <wuquanfeng@ihep.ac.cn>
# Feb 10, 2024
"""
    run_FORM(
        form_script_str::String;
        file_name::String="debug",
        multi_thread_flag::Bool=false
    )::String
  
Run the FORM script `form_script_str` and return the result.
"""
function run_FORM(
  form_script_str::String;
  file_name::String="debug",
  multi_thread_flag::Bool=false
)::String
  result_io = IOBuffer()

  try
      (run∘pipeline)(
          multi_thread_flag ?
              `$(tform()) -w$(Threads.nthreads()) -q -` :
              `$(form()) -q -`;
          stdin=IOBuffer(form_script_str),
          stdout=result_io
      ) # end run∘pipeline
  catch
    write( "$file_name.frm", form_script_str )
    @warn "Please check the $(joinpath(pwd(), "$file_name.frm"))!"
    rethrow()
  end # try

  return (String∘take!)( result_io )
end # function run_FORM










##############################################################################
# Changed by Quan-feng WU (wuquanfeng@ihep.ac.cn)
# March 9, 2023
"""
    make_amp_contraction_noexpand_script( 
        expr::Basic,
    )::String

Prepare the FORM script for the amplitude contraction, but do not do the expansion for the amplitude.
"""
function make_amp_contraction_noexpand_script(
    expr::Basic;
    dir::String=pwd()
)::String
##############################################################################

  model_parameters_content = (String∘read∘joinpath)( dir, "model_parameters.frm" )

  result_str = """
  #-

  Off Statistics;
  Off FinalStats;

  $(model_parameters_content)
  #include $(joinpath( art_dir(), "scripts", "contractor.frm" ))

  symbol sqrteta;

  format nospaces;
  format maple;

  Local expression = $(expr);
  .sort
  id GAij(spa1?,spa2?,mom?,mass?) = GAij(spa1,spa2,mom+mass*unity);
  .sort

  #call SimplificationNoExpand();

  #call contractDiracIndicesNoExpand();

  #call SimplificationNoExpand();

  #include $(joinpath( dir, "kin_relation.frm" ))
  .sort

  ***repeat;
  ***  id once FermionChain(?vars1, GA(mom?), ?vars2 ) = FV(mom,rho100)*FermionChain(?vars1, GA(rho100), ?vars2 );
  ***  sum rho100;
  ***endrepeat;


  id FV(mom?,rho?)*VecEpsilon?{VecEp,VecEpC}(int?,rho?,mom?,mass?) = 0;
  .sort

  while( match(FermionChain(?vars1,GA(rho?NonEPMU\$LORENTZ),?vars2)) );
    sum \$LORENTZ;
  endwhile;
  .sort
  *
  * Replace system dummy indices Nm_? by our dummy indices dum in case to read back to GiNaC.
  * We assume this should give the canonical form of FermionChain, 
  *   since it seems dummy indices Nm_? can make canonical form of an expression automatically.
  *

  ***repeat;
  ***if( match( SP(mom1?{q1,q2,q3,q4}\$MOM1,mom2?\$MOM2) ) );
  ***  id once SP(\$MOM1,\$MOM2) = FV(\$MOM1,rho1)*FV(\$MOM2,rho2)*LMT(rho1,rho2);
  ***  sum rho1;
  ***  sum rho2;
  ***endif;
  ***endrepeat; 
  ***.sort


  #do MUIDX = 1, 20, 1
    Multiply replace_(N`MUIDX'_?,dum`MUIDX');
  #enddo
  .sort

  ***id FV(rho1?,rho2?) = FV(rho1,rho2);
  ***id SP(rho1?,rho2?) = SP(rho1,rho2);
  ***.sort

  #write "%E", expression
  .sort

  .end

  """

  return result_str

end # function make_amp_contraction_noexpand_script










#############################################################
# Changed by Quan-feng WU (wuquanfeng@ihep.ac.cn)
# March 9, 2023
"""
    make_color_script( 
        color_factor::Basic,
    )::String

Specifically calculate the color factor `color_factor`.
"""
function make_color_script(
    color_factor::Basic
)::String
#############################################################

  result_str = """
  #-
  Off Statistics;
  
  format nospaces;
  format maple;
  
  #include $(joinpath( art_dir(), "scripts", "color.frm" ))
  
  Local colorFactor = $(color_factor);
  
  #call calc1_CF();
  .sort 
  
  #call calc2_CF();
  .sort 
  
  #write "%E", colorFactor
  .sort
  
  .end
  """

  return result_str

end # function make_color_script









