# Example julia -t 6 scripts/run_ardca.jl /Users/alex/proteins/aflatent/data/family/cm/cm_match_msa.fa testardca --verbose --lambdaJ 0.01 --lambdaH 0.01
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
            help = "reweighting threshold; defaults to auto (which corresponds to ??)"
            default = :auto
            arg_type = Union{Symbol, Float64}
        "--max_gap_fraction"
            help = "maximum fraction of gaps in a single sequence"
            arg_type = Float64
            default = 0.9
        "--permorder"  # TODO apply choices restriction, c.f. https://github.com/carlobaldassi/ArgParse.jl/issues/39
            help = "permutation of MSA columns"
            default = "entropic"
        "--lambdaJ"
            help = "regularisation strength for couplings"
            arg_type = Float64
            default = 0.02
        "--lambdaH"
            help = "regularisation strength for fields"
            arg_type = Float64
            default = 0.001
        "--maxit"
            help = "maximum number of iterations"
            arg_type = Int
            default = 1000
        "--verbose"
            help = "Run verbosely"
            action = :store_true
    end
    return parse_args(s)
end


parsed_args = parse_commandline()
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
)
