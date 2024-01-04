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

function check_form()::Tuple{Bool,String}
    if haskey( ENV, "FORM" ) && isfile( ENV["FORM"] )
        return true, ENV["FORM"]
    end # if
    return false, ""
end # function check_form

function check_tform()::Tuple{Bool,String}
    if haskey( ENV, "TFORM" ) && isfile( ENV["TFORM"] )
        return true, ENV["TFORM"]
    end # if
    return false, ""
end # function check_tform

function form()
    form_flag, form_path = check_form()
    form_flag && return `$form_path`
    return FORM_jll.form()
end

function tform()
    tform_flag, tform_path = check_tform()
    tform_flag && return `$tform_path`
    return FORM_jll.tform()
end

function build_qgraf(FC::String="gfortran")::String
    FC_flag, FC_path = check_FC(FC)
    @assert FC_flag "$FC not found."

    qgraf_installed_dir = (mkpath∘joinpath)( first(DEPOT_PATH), "local", "bin" )
    qgraf_installed_path = joinpath( qgraf_installed_dir, "qgraf" )
    
    qgraf_source = (first∘filter)( endswith(".f08"), readdir(Pkg.artifact"QGRAF"; join=true) )
    fmodules_dir = (mkdir∘joinpath)( mktempdir(), "fmodules" )

    main_file = (first ∘ filter)( startswith("main"), qgraf_source )
    license_header = readlines( joinpath( Pkg.dir("QGRAF"), main_file ) )[begin:93]
    printstyled( "Before building QGRAF, please read the following license:\n", color=:yellow)
    println( join(license_header, "\n") )
    printstyled( "If you agree with the license, please type \"yes\" to continue: ", color=:yellow )
    readline() ∈ ["yes", "YES", "Y", "y", "Yes"] || error("You must agree with the license to continue.")

    run(`$(FC_path) -o $(qgraf_installed_path) -Os -J $(fmodules_dir) $(qgraf_source)`)

    return joinpath( qgraf_installed_path )
end # function build_qgraf

function qgraf(;FC="gfortran")::Cmd
    qgraf_flag, qgraf_path = check_qgraf()
    qgraf_flag && return `$qgraf_path`
    return `$(build_qgraf(FC))`
end # function qgraf
