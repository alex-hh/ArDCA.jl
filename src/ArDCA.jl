module ArDCA

using ArgParse
using Random: randperm 
using Distributed: @distributed   
using Printf: @printf 
using LinearAlgebra: rmul!, norm, cholesky, inv
using ExtractMacro: @extract
using NLopt: Opt,ftol_abs!,xtol_rel!,xtol_abs!,ftol_rel!,maxeval!,min_objective!,optimize
using Distributions: wsample
using LoopVectorization: @avx 
using DCAUtils: read_fasta_alignment,remove_duplicate_sequences,compute_weights,compute_weighted_frequencies,add_pseudocount,compute_DI_gauss,compute_FN
using DCAUtils.ReadFastaAlignment: letter2num
using NPZ: npzread, npzwrite

# this determines which commands can be invoked from the global namespace in a session
# I guess non-exports can still be invoked via ArDCA.<name>??
export ardca,ArVar,ArAlg,ArNet,sample,sample_with_weights,sample_subsequence,epistatic_score,dms_single_site,loglikelihood,load_arnet,load_conf,write_vector,write_tuple_vector,gDCA,printrank,parse_commandline

include("types.jl")
include("ar.jl")
include("utils.jl")
include("dca.jl")
include("args.jl")


# entrypoint for app with packagecompiler
# https://julialang.github.io/PackageCompiler.jl/dev/apps.html#Creating-an-app
# https://github.com/JuliaLang/PackageCompiler.jl/tree/master/examples/MyApp
function julia_main()::Cint
    if isempty(ARGS)
        println("Must specify module name")
        return 1
    end

    try
        module_name=popfirst!(ARGS)
        if module_name == "gdca"
            # TODO specialise parse_commandline to gdca
            score_types = Dict("frob" => :frob, "DI" => :DI)

            parsed_args = parse_commandline()
            if parsed_args["autotheta"]
                parsed_args["theta"] = :auto
            end

            J, preds = gDCA(
                parsed_args["fastafile"],
                θ=parsed_args["theta"],
                max_gap_fraction=parsed_args["max_gap_fraction"],
                min_separation=parsed_args["min_separation"],
                score=score_types[parsed_args["score_type"]]
            )
            write_tuple_vector(preds, parsed_args["output_file"], colnames=("i", "j", "score"))
            if parsed_args["save_J"]
                println(size(J))
                contacts_file = string(split(parsed_args["output_file"],".")[1], ".npy")
                println(contacts_file)
                npzwrite(contacts_file, J)  # n.b this will be a bug if filename contains .
            end
            return 0

        else
            @show ARGS
            println("First argument must be module_name")
            println("Currently supported modules: gdca")
            return 1
        end

    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end

end


end # end module
