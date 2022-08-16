# Example julia -t 5 scripts/run_gdca.jl data/PF14/PF00014_mgap6.fasta outputs/ardca_PF14_conreg --verbose
using ArgParse
using ArDCA
using NPZ

score_types = Dict("frob" => :frob, "DI" => :DI)

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "fastafile"
            help = "Input MSA"
            required = true
        "output_file"
            help = "File in which to store predicted contacts"
            required = true
        "--theta"
            help = "1 - reweighting pairwise identity threshold"
            default = 0.2
            arg_type = Float64
        "--autotheta"
            help = "automatically compute reweighting threshold. Computed as mean value of the similarity fraction between all possible pairs of sequences by dependency DCAUtils"
            action = :store_true
        "--max_gap_fraction"
            help = "maximum fraction of gaps in a single sequence"
            arg_type = Float64
            default = 0.9
        "--min_separation"
            help = "Minimum contact separation"
            arg_type = Int
            default = 1
        "--save_J"
            help="save inverse covariance matrix as np array"
            action = :store_true
        "--pseudocount"
            help="pseudocount (0,1). 0.8 recommended if using frob scoring, 0.2 for DI"
            default = 0.8
            arg_type = Float64
        "--score_type"
            help="method for computing predicted contacts from inverse covariance (frob or DI)"
            default = "frob"
            arg_type = String
    end
    return parse_args(s)
end


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
