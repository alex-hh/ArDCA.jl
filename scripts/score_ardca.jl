# score a sequence input directly, printing the result to stdout
# realistically we want to be able to score a batch of sequences at once
# the best way of doing this is probably (unfortunately) going to be
# to write the scoring logic directly in python. (or poss using PyCall/PyJulia)
# https://discourse.julialang.org/t/is-pyjulia-being-maintained/82607/40
using ArDCA
using ArgParse
using DCAUtils


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "ardca_filename"
            help = "Trained ArDCA parameters/conf base name"
            required = true
        "--sequence"
        	help = "Sequence to score"
        	arg_type = String
        "--fasta_file"
            help = "Fasta to score"
            arg_type = String
        "--output_file"
            help = "Text file in which to save scores"
            arg_type = String
    end
    return parse_args(s)
end


# TODO check sequence is correct length somehow?
parsed_args = parse_commandline()
arnet = load_arnet(parsed_args["ardca_filename"])
if !isnothing(parsed_args["sequence"])
    val = loglikelihood(parsed_args["sequence"], arnet)
    println(val)
    if !isnothing(parsed_args, ["output_file"])
        io = open(filename, "w") do io
            println(io, val)
        end
    end
elseif !isnothing(parsed_args["fasta_file"])
    Z = read_fasta_alignment(parsed_args["fasta_file"],  1.0)  # second arg is max fraction of gaps
    vals = loglikelihood(Z, arnet)
    if !isnothing(parsed_args["output_file"])
        write_vector(vals, parsed_args["output_file"])
    end
else
    error("Must pass either an input sequence or fasta file to score")
end
