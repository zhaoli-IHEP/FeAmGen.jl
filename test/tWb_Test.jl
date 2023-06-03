using Dates, FeAmGen, FeynUtils 

start = now()
@info "tWb_Test starts @ $(start)"

#----------------------------------------------------------------------------
# t->b+W 1-loop, 2-loop including tbW vertices tests
#----------------------------------------------------------------------------
seed_proc_yaml_str( ; nloop::Int64 = 2::Int64 ) = """
# input file for calculation details
# model related information

# model name
model_name: "sm_CKMdiag_Haa"

# use unitary_gauge for internal vector boson
unitary_gauge: true

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
Amp_QCD_order: $(2*nloop) 
# order of QED coupling ee in the amplitude
Amp_QED_order: 1  
# order of special coupling in the amplitude
Amp_SPC_order: 0  

# min ep power in the amplitude
Amp_Min_Ep_Xpt: $(-2*nloop)
# max ep power in the amplitude
Amp_Max_Ep_Xpt: 0

# incoming and outgoing information
incoming: [ "t" ]          # incoming particles
outgoing: [ "Wplus", "b" ]               # outgoing particles 

# Symmetry configuration
symmetry: []

"""


#-------------------------------------------
# Start running
for nloop in [3]

  open( "seed_tWb_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, seed_proc_yaml_str(nloop=nloop) )
  end # close

  card_list = digest_seed_proc( "seed_tWb_proc_$(nloop)Loop.yaml" )

  for one_card in card_list
    box_message( "[ Generate amplitudes for $(one_card) ]" )
    generate_amp( one_card )
  end # for one_card

end # for nloop

@info "tWb_Test ends @ $(now()) started from $(start)"


