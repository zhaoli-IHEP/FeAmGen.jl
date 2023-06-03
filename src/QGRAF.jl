# check if qgraf
function check_qgraf()::Tuple{Bool,String}
    if haskey( ENV, "QGRAF" ) && isfile( ENV["QGRAF"] )
        return true, ENV["QGRAF"]
    end # if
    io = IOBuffer()
    p = (run∘pipeline)( Cmd( `which qgraf`; ignorestatus=true ), stdout=io )
    iszero( p.exitcode ) && return true, (String∘take!)(io)[begin:end-1]

    builded_qgraf = joinpath( first(DEPOT_PATH), "local", "bin", "qgraf" )
    isfile( builded_qgraf ) && return true, builded_qgraf
    
    return false, ""
end # function check_qgraf

function check_FC(FC::String="gfortran")::Tuple{Bool,String}
    if haskey( ENV, "FC" ) && isfile( ENV["FC"] )
        return true, ENV["FC"]
    end # if
    io = IOBuffer()
    p = (run∘pipeline)( Cmd( `which $FC`; ignorestatus=true ), stdout=io )
    iszero( p.exitcode ) && return true, (String∘take!)(io)[begin:end-1]
    
    return false, ""
end # function check_FC

function build_qgraf(FC::String="gfortran")::String
    FC_flag, FC_path = check_FC(FC)
    @assert FC_flag "$FC not found."

    qgraf_installed_dir = (mkpath∘joinpath)( first(DEPOT_PATH), "local", "bin" )
    qgraf_installed_path = joinpath( qgraf_installed_dir, "qgraf" )
    
    qgraf_source = (first∘filter)( endswith(".f08"), readdir(Pkg.artifact"QGRAF"; join=true) )
    fmodules_dir = (mkdir∘joinpath)( mktempdir(), "fmodules" )
    run(`$(FC_path) -o $(qgraf_installed_path) -Os -J $(fmodules_dir) $(qgraf_source)`)

    return joinpath( qgraf_installed_path )
end # function build_qgraf

function qgraf(;FC="gfortran")::Cmd
    qgraf_flag, qgraf_path = check_qgraf()
    qgraf_flag && return `$qgraf_path`
    return `$(build_qgraf(FC))`
end # function qgraf
