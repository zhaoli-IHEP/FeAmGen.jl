using Dates, FeAmGen

@info "ggttbar_CT_Test starts @ $(now())"


#----------------------------------------------------------------------------
# gg->ttbar with counter terms
#----------------------------------------------------------------------------
generic_ggttbar_seed_proc_yaml_str( ; n_loop::Int64=2, CT_order::Int64=0 ) = """
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
n_loop: $(n_loop - CT_order)
# order of QCD counter-term vertices
QCDCT_order: $(CT_order)

# order of QCD coupling gs in the amplitude
Amp_QCD_order: $(2 + 2*n_loop) 
# order of QED coupling ee in the amplitude
Amp_QED_order: 0  
# order of special coupling in the amplitude
Amp_SPC_order: 0  

# min ep power in the amplitude
Amp_Min_Ep_Xpt: $(-2*n_loop)
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
n_loop = 2
CT_order = 1

write( 
  "seed_ggttbar_proc_$(n_loop)Loop_with_$(CT_order)CT.yaml",
  generic_ggttbar_seed_proc_yaml_str(n_loop=n_loop, CT_order=CT_order)
) # write

card_list = digest_seed_proc( "seed_ggttbar_proc_$(n_loop)Loop_with_$(CT_order)CT.yaml" )
filter!( contains("g_g_TO_t_tbar.yaml"), card_list )
@assert !isempty(card_list)
  
# we only test the gluon fusion subprocess
generate_amp.( card_list )

@info "ggttbar_Test ends @ $(now())"


