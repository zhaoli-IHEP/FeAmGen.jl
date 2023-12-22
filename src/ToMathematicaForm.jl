# Copyright (c) 2023 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

function to_m_file(amp_file::String, m_file::Union{Missing, String}=missing)::Nothing
    @assert isfile(amp_file) "Please check the file path of $amp_file!"
    amp_jld = jldopen(amp_file, "r")
    mma_file = !ismissing(m_file) ? m_file : begin
        filename, _ = splitext(amp_file)
        filename * ".m"
    end

    loop_den_list = to_Mathematica_form(amp_jld["loop_den_list"])
    loop_den_xpt_list = amp_jld["loop_den_xpt_list"]

    kin_relation = Dict{String, String}()
    for (key, value) ∈ amp_jld["kin_relation"]
        key_mma_form = to_Mathematica_form(key)
        value_mma_form = to_Mathematica_form(value)
        kin_relation[key_mma_form] = value_mma_form
    end

    signed_symmetry_factor = amp_jld["signed_symmetry_factor"]

    amp_color_list = to_Mathematica_form(amp_jld["amp_color_list"])
    amp_lorentz_list = to_Mathematica_form(amp_jld["amp_lorentz_list"])

    mom_symmetry = amp_jld["mom_symmetry"]
    color_symmetry = amp_jld["color_symmetry"]

    model_coupling_dict = Dict{String, String}()
    for (key, value) ∈ amp_jld["model_coupling_dict"]
        key_mma_form = to_Mathematica_form(key)
        value_mma_form = to_Mathematica_form(value)
        model_coupling_dict[key_mma_form] = value_mma_form
    end

    model_parameter_dict = Dict{String, String}()
    for (key, value) ∈ amp_jld["model_parameter_dict"]
        key_mma_form = to_Mathematica_form(key)
        value_mma_form = to_Mathematica_form(value)
        model_parameter_dict[key_mma_form] = value_mma_form
    end

    open(mma_file, "w+") do io
        write(io, "ampColorList = {" * join(amp_color_list, ", ") * "};\n")
        write(io, "ampLorentzList = {" * join(amp_lorentz_list, ", ") * "};\n")
        write(io, "\n")
        write(io, "loopDenList = {" * join(loop_den_list, ", ") * "};\n")
        write(io, "loopDenXptList = {" * join(loop_den_xpt_list, ", ") * "};\n")
        write(io, "\n")
        write(io, "signedSymmetryFactor = $signed_symmetry_factor;\n")
        write(io, "\n")
        write(io,
            "kinRelation = {" * join(
                ["$key -> $value" for (key, value) ∈ kin_relation], ", "
            ) * "};\n"
        )
        write(io, "\n")
        write(io,
            "momSymmetry = {" * join(
                ["$key -> $value" for (key, value) ∈ mom_symmetry], ", "
            ) * "};\n"
        )
        write(io, "colorSymmetry = {" * join(
                ["$key -> $value" for (key, value) ∈ color_symmetry], ", "
            ) * "};\n"
        )
        write(io, "\n")
        write(io,
            "modelCouplingDict = {" * join(
                ["$key -> $value" for (key, value) ∈ model_coupling_dict], ", "
            ) * "};\n"
        )
        write(io,
            "modelParameterDict = {" * join(
                ["$key -> $value" for (key, value) ∈ model_parameter_dict], ", "
            ) * "};\n"
        )
    end

    return nothing
end

function run_FORM(
    form_script_str::String;
    file_name::String="debug",
    multi_thread_flag::Bool=false
)::String
    result_io = IOBuffer()

    try
        (run ∘ pipeline)(
            multi_thread_flag ?
                `$(tform()) -w$(Threads.nthreads()) -q -` :
                `$(form()) -q -`;
            stdin=IOBuffer(form_script_str),
            stdout=result_io
        )
    catch
        write("$file_name.frm", form_script_str)
        @warn "Please check the $(joinpath(pwd(), "$file_name.frm"))!"
        rethrow()
    end

    return (String ∘ take!)(result_io)
end

to_Mathematica_form(input::Basic)::String = (to_Mathematica_form ∘ string)(input)
function to_Mathematica_form(input_str::String)::String
    
    floating_str_list, input = begin
        floating_pattern = r"(\d+)(\.\d+)"
        floating_range_list = findall(floating_pattern, input_str)
        floating_str_list = [input_str[floating_range] for floating_range ∈ floating_range_list]
        for (floating_index, floating_str) ∈ enumerate(floating_str_list)
            input_str = replace(input_str, floating_str => "floatingSymbol$floating_index"; count=1)
        end
        floating_str_list, Basic(input_str)
    end

    old_free_symbol_names = map(string, free_symbols(input))
    old_function_names = map(get_name, function_symbols(input))
    unique!(old_function_names)

    new_free_symbol_names = String[]
    new_function_names = String[]
    input_str = string(input)
    for old_name ∈ old_free_symbol_names
        new_name = replace(old_name, "_" => "")
        duplicate_time = 0
        original_new_name = new_name
        while new_name ∈ new_free_symbol_names
            duplicate_time += 1
            new_name = original_new_name * string(duplicate_time)
            @warn """
            Duplicate symbol name $original_new_name from $(old_name)!
            We will add a index to the new symbol name as $new_name.
            """
        end
        push!(new_free_symbol_names, new_name)
        input_str = replace(input_str, old_name => new_name)
    end
    for old_name ∈ old_function_names
        new_name = replace(old_name, "_" => "")
        duplicate_time = 0
        original_new_name = new_name
        while new_name ∈ new_function_names
            duplicate_time += 1
            new_name = original_new_name * string(duplicate_time)
            @warn """
            Duplicate function name $original_new_name from $(old_name)!
            We will add a index to the new function name as $new_name.
            """
        end
        push!(new_function_names, new_name)
        input_str = replace(input_str, old_name => new_name)
    end

    isempty(new_free_symbol_names) && push!(new_free_symbol_names, "noSymbol")
    isempty(new_function_names) && push!(new_function_names, "noFunction")

    Form_script_str = """
    #-

    Format Mathematica;
    Format nospaces;

    Off Statistics;
    Off FinalStats;

    Symbol pi, $(join(new_free_symbol_names, ", "));
    CFunction sqrt, $(join(new_function_names, ", "));

    Local expr = $input_str;
    .sort

    #write "%E", expr
    .sort

    .end
    """

    result_str = run_FORM(Form_script_str)
    result_str = replace(result_str,
        " " => "",
        '\n' => "",
        "conj[" => "Conjugate[",
        "sqrt[" => "Sqrt[",
        "im" => "I",
        "pi" => "Pi"
    )
    result_str = replace(result_str,
        [
            "floatingSymbol$floating_index" => floating_str
                for (floating_index, floating_str) ∈ enumerate(floating_str_list)
        ]...
    )

    return result_str
end
to_Mathematica_form(input::Array) = map(to_Mathematica_form, input)
