using Dates, FeAmGen, FeynUtils 

@info "ggttbar_Test starts @ $(now())"


#----------------------------------------------------------------------------
# gg->ttbar 0-loop, 1-loop tests
#----------------------------------------------------------------------------
generic_ggttbar_seed_proc_yaml_str( ; nloop::Int64 = 2::Int64 ) = """
# input file for calculation details
# model related information

# model name
model_name: "sm"

# use unitary_gauge for internal vector boson
unitary_gauge: false

# content of "quark-parton", the "parton" will also contain gluon. Only for seed program.
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "s", "c", "b", "sbar", "cbar", "bbar" ] )   
# for single top @ NNLO
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "b", "bbar" ] )   
# for Higgs+Jet @ NLO
partons: [ "g", "u", "ubar", "d", "dbar", "b", "bbar" ] 

# only for seed program
AllowLeptonNumberViolation: false
AllowQuarkGenerationViolation: false

# process information
DropTadpole: true              # drop tadpole?
DropWFcorrection: true         # drop WFcorrection?

# number of loops
n_loop: $(nloop)
# order of QCD counter-term vertices
QCDCT_order: 0   

# order of QCD coupling gs in the amplitude
Amp_QCD_order: $(2+2*nloop) 
# order of QED coupling ee in the amplitude
Amp_QED_order: 0  
# order of special coupling in the amplitude
Amp_SPC_order: 0  

# min ep power in the amplitude
Amp_Min_Ep_Xpt: $(-2*nloop)
# max ep power in the amplitude
Amp_Max_Ep_Xpt: 0

# incoming and outgoing information
incoming: [ "parton", "parton" ]          # incoming particles
outgoing: [ "t", "tbar" ]               # outgoing particles 

# Symmetry configuration
momentum_symmetry: []
color_symmetry: []

"""


#-------------------------------------------
# Start running
for nloop in [0,1]

  open( "seed_ggttbar_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, generic_ggttbar_seed_proc_yaml_str(nloop=nloop) )
  end # close

  card_list = digest_seed_proc( "seed_ggttbar_proc_$(nloop)Loop.yaml" )
  @assert "parton_parton_TO_t_tbar_$(nloop)Loop/g_g_TO_t_tbar.yaml" in card_list
  
  # we only test the gluon fusion subprocess
  generate_amp( "parton_parton_TO_t_tbar_$(nloop)Loop/g_g_TO_t_tbar.yaml" )

end # for nloop


@info "ggttbar_Test ends @ $(now())"


