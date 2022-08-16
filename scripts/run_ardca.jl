# Example julia -t 5 scripts/run_gdca.jl data/PF14/PF00014_mgap6.fasta outputs/ardca_PF14_conreg --verbose
using ArDCA
using ArgParse


permorders = Dict("natural" => :NATURAL, "entropic" => :ENTROPIC, "rev_entropic" => :REV_ENTROPIC, "random" => :RANDOM)


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "fastafile"
            help = "Input MSA"
            required = true
        "output_file"
        	help = "File in which to store trained parameters"
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
        "--permorder"  # TODO apply choices restriction, c.f. https://github.com/carlobaldassi/ArgParse.jl/issues/39
            help = "permutation of MSA columns"
            default = "entropic"
        "--lambdaJ"
            help = "regularisation strength for couplings (paper uses 0.01 for contact experiments, 0.0001 for generation)"
            arg_type = Float64
            default = 0.01
        "--lambdaH"
            help = "regularisation strength for fields (paper uses 0.0001 for contact experiments, 0.000001 for generation)"
            arg_type = Float64
            default = 0.0001
        "--maxit"
            help = "maximum number of iterations"
            arg_type = Int
            default = 1000
        "--verbose"
            help = "Run verbosely"
            action = :store_true
        "--save_fp16"
            help = "Save params at half precision"
            action = :store_true
        "--save_txt"
            help = "Save text file containing params"
            action = :store_true
    end
    return parse_args(s)
end


parsed_args = parse_commandline()
if parsed_args["autotheta"]
    parsed_args["theta"] = :auto
end

arnet,arvar=ardca(
	parsed_args["fastafile"],
	theta=parsed_args["theta"],
	max_gap_fraction=parsed_args["max_gap_fraction"],
	verbose=parsed_args["verbose"],
	lambdaJ=parsed_args["lambdaJ"],
	lambdaH=parsed_args["lambdaH"],
	output_file=parsed_args["output_file"],
	maxit=parsed_args["maxit"],
	permorder=permorders[parsed_args["permorder"]],
    save_fp16=parsed_args["save_fp16"]
)
