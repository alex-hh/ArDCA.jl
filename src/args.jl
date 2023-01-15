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