# Draw samples from a trained network and optionally compute statistics.
using ArDCA
using ArgParse
using DCAUtils


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "ardca_filename"
            help = "Trained ArDCA parameters/conf base name"
            required = true
        "output_file"
            help = "File in which to save samples"
            required = true
        "--reference_sequence"
            help = "Reference sequence, if not specified use sequence in training msa (index specified by train_sequence_index)"
            default = nothing
            arg_type = Union{String, Nothing}
        "--msa"
            help = "Training MSA"
            default = nothing
            arg_type = Union{String, Nothing}
        "--train_sequence_index"
            help = "Index of sequence in training MSA to use as reference to compute epistatic score"
            arg_type = Int
            default = 1
        "--min_separation"
            help = "Minimum sequence separation between contacts"
            arg_type = Int
            default = 1
    end
    return parse_args(s)
end


parsed_args = parse_commandline()
arnet = load_arnet(parsed_args["ardca_filename"])
cfg_dict = load_conf(string(parsed_args["ardca_filename"], ".conf"))
if !isnothing(parsed_args["reference_sequence"])
    contact_scores = epistatic_score(
        arnet,
        parsed_args["reference_sequence"],
        min_separation=parsed_args["min_separation"],
    )
elseif !isnothing(parsed_args["msa"])
    max_gap_fraction = 0.5
    Z = read_fasta_alignment(parsed_args["msa"], max_gap_fraction)
    contact_scores = epistatic_score(
        arnet,
        Z,
        parsed_args["train_sequence_index"],
        min_separation=parsed_args["min_separation"],
    )
else
    throw("must pass either reference sequence or MSA")
end

# contact scores is a tuple: (i, j, frobnorm), we want to write this to csv.
write_tuple_vector(contact_scores, parsed_args["output_file"], colnames=("i", "j", "score"))
