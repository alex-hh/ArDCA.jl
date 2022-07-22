# Draw samples from a trained network and optionally compute statistics.
using ArDCA
using ArgParse


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "ardca_filename"
            help = "Trained ArDCA parameters/conf base name"
            required = true
        "output_file"
        	help = "File in which to save samples"
        	required = true
        "--num_sequences"
        	help = "Number of sequences to sample"
        	arg_type = Int
        	default = 10000
        "--compute_weights"
        	help = "Compute log likelihood of samples"
        	action = :store_true
        "--seed_sequence"
        	default = nothing
        	arg_type = Union{String, Nothing}
        	help = "Seed sequence"
        "--no_gaps"
        	help = "If this flag is passed, we won't sample gaps"
        	action = :store_true
        "--temperature"
        	help = "Softmax temperature"
        	arg_type = Float64
        	default = 1.0
    end
    return parse_args(s)
end


# https://github.com/pagnani/ArDCA.jl/issues/19
function write_fasta(filedest::String, Z; W::Union{Vector{Float64},Nothing}=nothing)
	num2let = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I','K', 'L','M', 'N','P','Q', 'R', 'S','T', 'V', 'W', 'Y','-']
	N,M = size(Z)
	open(filedest,"w") do fp
   		for s in 1:M
       	  	if isnothing(W)
       	    	label = ">Seq$s"
       	    else
       	    	label = string(">Seq$s", "_", "logp", log(W[s]))
       	    end
			println(fp, label)
			for a in 1:N
            	print(fp,num2let[Int(Z[a,s])])
			end
				println(fp)
     	end
   	end
end


parsed_args = parse_commandline()
arnet = load_arnet(parsed_args["ardca_filename"])

W = nothing
seed = parsed_args["seed_sequence"]
M = parsed_args["num_sequences"]
tau = parsed_args["temperature"]
if !isnothing(seed)
	# TODO check how loglikelihood is handled in this case: does it include seed? A. I believe yes.
	# TODO check how seed works - does it assume 'natural' ordering (for seed itself?)
	W, gen = sample_subsequence(seed, arnet, M, no_gaps=parsed_args["no_gaps"], temperature=tau)
elseif parsed_args["compute_weights"]
	W, gen = sample_with_weights(arnet, M, no_gaps=parsed_args["no_gaps"], temperature=tau)
else
	gen = sample(arnet, M, no_gaps=parsed_args["no_gaps"], temperature=tau)
end

write_fasta(parsed_args["output_file"], gen, W=W)
