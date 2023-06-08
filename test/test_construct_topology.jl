using FeAmGen

println()
println( "--"^14 )
@info "Test for `construct_den_topology`"

rootdir = dirname(@__FILE__)
#proc_dir = joinpath( rootdir, "t_TO_t_4Loop" )
#amp_dir = joinpath( proc_dir, "t_TO_t_4Loop_amplitudes" )
#proc_dir = joinpath( rootdir, "Wplus_t_TO_Wplus_t_3Loop" )
#amp_dir = joinpath( proc_dir, "Wplus_t_TO_Wplus_t_3Loop_amplitudes" )
proc_dir = joinpath( rootdir, "b_g_TO_Wminus_t_2Loop" )
amp_dir = joinpath( proc_dir, "b_g_TO_Wminus_t_2Loop_amplitudes" )

@assert isdir(amp_dir)
FeAmGen.construct_den_topology( amp_dir )


