using FeAmGen

println()
println( "--"^14 )
@info "Test for `construct_den_topology`"

proc_str_list = [ "b_g_TO_Wminus_t_2Loop", "Wplus_t_TO_Wplus_t_3Loop", "t_TO_t_4Loop" ]
proc_str = proc_str_list[3]

rootdir = dirname(@__FILE__)
proc_dir = joinpath( rootdir, proc_str )
amp_dir = joinpath( proc_dir, "$(proc_str)_amplitudes" )

@assert isdir(amp_dir)
construct_den_topology( amp_dir, mom_shift_opt = true )


