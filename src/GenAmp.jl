
######################################################################
"""
    generate_amp( proc_file::String )::Nothing

Read-in the card for the specific process and produce the relevant amplitude.
This function is one of the front-end functions in FeAmGen.jl.
The directory of model files are supposed in ".".
"""
function generate_amp( proc_file::String; model_paths=[pwd()] )::Nothing
######################################################################

  #------------------------------------------------------------------
  @assert isfile(proc_file) "The first argument is not a file!"
  working_dir = pwd()
  proc_dir = begin
    tmp_dir = dirname(proc_file)
    isempty(tmp_dir) ? pwd() : tmp_dir
  end # proc_dir
  input = YAML.load_file( proc_file )
  #------------------------------------------------------------------


  #------------------------------------------------------------------
  @info "Choose model" model=input["model_name"]
  model = readin_model( input; model_paths=model_paths )
  cd(proc_dir)
  generate_QGRAF_model( model )
  logging_model( model )
  #------------------------------------------------------------------

  @info "Usage of unitary gauge" unitary_gauge=input["unitary_gauge"]

  @info "Drop Tadpole" DropTadpole=input["DropTadpole"]

  @info "Drop WFcorrection" DropWFcorrection=input["DropWFcorrection"]

  @info "Number of loops" n_loop=input["n_loop"]

  @info "QCD CT-order" QCDCT_order=input["QCDCT_order"]

  @info "Order of QCD coupling gs in the amplitude" Amp_QCD_order=input["Amp_QCD_order"]

  @info "Order of QED coupling ee in the amplitude" Amp_QED_order=input["Amp_QED_order"]

  @info "Order of special coupling in the amplitude" Amp_SPC_order=input["Amp_SPC_order"]

  @info "Min ep power in the amplitude" Amp_Min_Ep_Xpt=input["Amp_Min_Ep_Xpt"]
  @info "Max ep power in the amplitude" Amp_Max_Ep_Xpt=input["Amp_Max_Ep_Xpt"]

  @info "Incoming" input["incoming"]
  @info "Outgoing" input["outgoing"]

  @info "Coupling factor" couplingfactor=input["couplingfactor"]

  #----------------------------------------------------------------------
  # Run the QGRAF
  generate_Feynman_diagram( model, input )

  generate_amplitude( model, input )

  cd(working_dir)

  return nothing

end # function generate_amp




