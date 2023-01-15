using ArgParse
using ArDCA
using NPZ


score_types = Dict("frob" => :frob, "DI" => :DI)


J, preds = gDCA(
    "data/PF14/PF00014_mgap6.fasta",
    θ=0.2,
    max_gap_fraction=0.9,
    min_separation=1,
    score=score_types["frob"],
)
write_tuple_vector(preds, "outputs/gdca_PF14_conreg.csv", colnames=("i", "j", "score"))
println(size(J))
params_file = string(split("outputs/gdca_PF14_conreg",".")[1], ".npy")
println(params_file)
npzwrite(params_file, J)  # n.b this will be a bug if filename contains .
