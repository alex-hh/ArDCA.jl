# Example julia -t 5 scripts/run_gdca.jl data/PF14/PF00014_mgap6.fasta outputs/ardca_PF14_conreg --verbose
using ArgParse
using ArDCA
using NPZ

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
