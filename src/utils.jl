# https://discourse.julialang.org/t/how-to-export-a-column-vector-to-text-file-with-writedlm-function/37818
function write_vector(vec::AbstractVector, filename::AbstractString)
    io = open(filename, "w") do io
        for x in vec
            println(io, x)
        end
    end
end

function read_vector(filename::AbstractString)
    vec = Vector{Float64}()
    io = open(filename, "r") do io
        for line in readlines(io)
            if !isempty(line)
                # println(line)
                append!(vec, parse(Float64, line))
            end
        end
    end
    return vec
end

function write_dict(d::Dict, filename::AbstractString)
    io = open(filename, "w") do io
        for (key, val) in d
            println(io, string(key, ",", val))
        end
    end
end

function read_dict(filename::AbstractString)
    d = Dict()
    io = open(filename, "r") do io
        for line in readlines(io)
            if !isempty(line)
                key, val = split(line, ",")
                d[key] = val
            end
        end
    end
    return d
end

function read_fasta(filename::AbstractString, max_gap_fraction::Real, theta::Any, remove_dups::Bool)
    Z = read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups
        Z, _ = remove_duplicate_sequences(Z)
    end
    N, M = size(Z)
    q = round(Int, maximum(Z))
    W, Meff = compute_weights(Z, theta)
    return W, Z, N, M, q
end

# these are empirical frequencies in the first column
function computep0(W, Z, q)
    # @extract var:W Z q
    p0 = zeros(q)
    for i in 1:length(W)
        p0[Z[1, i]] += W[i]
    end
    p0
end

function compute_empirical_freqs(Z::AbstractArray{Ti,2}, W::AbstractVector{Float64}, q::Ti) where {Ti<:Integer}
    N, M = size(Z)
    f = zeros(q, N)
    @inbounds for i in 1:N
        for s in 1:M
            f[Z[i, s], i] += W[s]
        end
    end
    f
end

function entropy(Z::AbstractArray{Ti,2}, W::AbstractVector{Float64}) where {Ti<:Integer}
    N, M = size(Z)
    q = maximum(Z)
    f = compute_empirical_freqs(Z, W, q)
    S = zeros(N)
    @inbounds for i in 1:N
        _s = 0.0
        for a in 1:q
            _s -= f[a, i] > 0 ? f[a, i] * log(f[a, i]) : 0.0
        end
        S[i] = _s
    end
    S
end

function load_conf(cfg_filename::String)
    permorders = Dict("NATURAL" => :NATURAL, "ENTROPIC" => :ENTROPIC, "REV_ENTROPIC" => :REV_ENTROPIC, "RANDOMus" => :RANDOM)
    d = read_dict(cfg_filename)
    d["lambdaJ"] = parse(Float64, d["lambdaJ"])
    d["lambdaH"] = parse(Float64, d["lambdaH"])
    d["permorder"] = permorders[d["permorder"]]  # convert to symbol
    d["epsconv"] = parse(Float64, d["epsconv"])
    d["maxit"] = parse(Int, d["maxit"])
    idxperm = map(x->parse(Int, x), split(d["idxperm"], ";"))
    p0 = map(x->parse(Float64,x), split(d["p0"], ";"))
    d["p0"] = p0
    d["idxperm"] = idxperm
    d["N"] = parse(Int, d["N"])
    d["q"] = parse(Int, d["q"])
    return d
end

function load_arnet(filename::String)
    params = read_vector(string(filename, ".params"))
    cfg_dict = load_conf(string(filename, ".conf"))
    arvar = ArPredVar(cfg_dict["N"], cfg_dict["q"], cfg_dict["idxperm"], cfg_dict["p0"])
    return ArNet(arvar.idxperm, unpack_params(params, arvar)...)  # ... here unpacks the tuple
end

function unpack_params(θ, arvar::Union{ArVar,ArPredVar})
    @extract arvar:N q p0

    arrJ = Array{Float64,3}[]
    arrH = Vector{Float64}[]
    ctr = 0
    for site in 1:N-1
        _arrJ = zeros(q, q, length(1:site))
        for i in 1:site
            for a in 1:q
                for b in 1:q
                    ctr += 1
                    _arrJ[b, a, i] = θ[ctr]
                end
            end
        end
        push!(arrJ, _arrJ)
        _arrH = zeros(q)
        for a in 1:q
            ctr += 1
            _arrH[a] = θ[ctr]
        end
        push!(arrH, _arrH)
    end
    @assert ctr == length(θ)
    p0, arrJ, arrH
end

function softmax!(r::Vector{Float64}, x::Vector{Float64})
    n = length(x)
    u = maximum(x)
    s = 0.0
    @inbounds @avx for i = 1:n
        r[i] = exp(x[i] - u)
        s += r[i]
    end
    invs = inv(s)
    @inbounds @avx for i = 1:n
        r[i] *= invs
    end
    r
end

"""
    softmax(x::AbstractArray{<:Real})
Return the [`softmax transformation`](https://en.wikipedia.org/wiki/Softmax_function) applied to `x`
"""
softmax!(x::Vector{Float64}) = softmax!(x, x)
softmax(x::Vector{Float64}) = softmax!(similar(x, Float64), x)

"""
    sample(arnet::ArNet, msamples::Int)
Return a generated alignment in the form of a `N × msamples`  matrix of type `::Matrix{Int}`  
# Examples
```
julia> arnet,arvar=ardca("file.fasta",verbose=true,permorder=:ENTROPIC, lambdaJ=0.001,lambdaH=0.001);
julia> Zgen=Zgen=sample(arnet,1000);
```
"""
function sample(arnet::ArNet, msamples::Int)
    @extract arnet:H J p0 idxperm
    q = length(p0)
    N = length(H) # here N is N-1 !!
    backorder = sortperm(idxperm)
    res = Matrix{Int}(undef, N + 1, msamples)
    Threads.@threads for i in 1:msamples
        totH = Vector{Float64}(undef, q)
        sample_z = Vector{Int}(undef, N + 1)
        sample_z[1] = wsample(1:q, p0)
        for site in 1:N
            Js = J[site]
            h = H[site]
            copy!(totH, h)
            @avx for i in 1:site
                for a in 1:q
                    totH[a] += Js[a, sample_z[i], i]
                end
            end
            p = softmax(totH)
            sample_z[site+1] = wsample(1:q, p)
        end
        res[:, i] .= sample_z
    end
    permuterow!(res, backorder)
end

"""
    sample_with_weights(arnet::ArNet, msamples::Int)
Return a generated alignment in the form of a `N × msamples`  matrix of type `::Matrix{Int}` and the relative probabilities under the module

# Examples
```
julia> arnet,arvar=ardca("file.fasta",verbose=true,permorder=:ENTROPIC, lambdaJ=0.001,lambdaH=0.001);
julia> Wgen,Zgen=sample_with_weights(arnet,1000);
```
"""
function sample_with_weights(arnet::ArNet, msamples)
    @extract arnet:H J p0 idxperm
    q = length(p0)
    N = length(H) # here N is N-1 !!
    backorder = sortperm(idxperm)
    W = Vector{Float64}(undef, msamples)
    res = Matrix{Int}(undef, N + 1, msamples)
    Threads.@threads for i in 1:msamples
        totH = Vector{Float64}(undef, q)
        sample_z = Vector{Int}(undef, N + 1)
        sample_z[1] = wsample(1:q, p0)
        logw = log(p0[sample_z[1]])
        for site in 1:N
            Js = J[site]
            h = H[site]
            copy!(totH, h)
            @avx for i in 1:site
                for a in 1:q
                    totH[a] += Js[a, sample_z[i], i]
                end
            end
            p = softmax(totH)
            sample_z[site+1] = wsample(1:q, p)
            logw += log(p[sample_z[site+1]])
        end
        W[i] = exp(logw)
        res[:, i] .= sample_z
    end
    return W, permuterow!(res, backorder)
end

"""
    sample_subsequence(x::String, arnet::ArNet, msamples)
Return a generated alignment in the form of a `N × msamples`  matrix of type
`::Matrix{Int}` and the relative probabilities under the model. The alignment is 
forced to start with with a sequence `x` (in amino acid single letter alphabet) 
and then autoregressively generated.

# Example
```
julia> arnet,arvar=ardca("file.fasta",verbose=true,permorder=:ENTROPIC, lambdaJ=0.001,lambdaH=0.001);
julia> Wgen,Zgen=sample_subsequence("MAKG",arnet,1000);
```
"""
function sample_subsequence(x::String, arnet::ArNet, msamples)
    x0 = letter2num.(collect(x))
    return sample_subsequence(x0, arnet, msamples)
end

"""
    sample_subsequence(x::Vector{T}, arnet::ArNet, msamples)
Return a generated alignment in the form of a `N × msamples`  matrix of type
`::Matrix{Int}` and the relative probabilities under the model. The alignment is 
forced to start with with a sequence `x` (in integer number coding) 
and then autoregressively generated.

# Example
```
julia> arnet,arvar=ardca("file.fasta",verbose=true,permorder=:ENTROPIC, lambdaJ=0.001,lambdaH=0.001);
julia> Wgen,Zgen=sample_subsequence([11,1,9,6],arnet,1000);
```
"""
function sample_subsequence(x0::Vector{T}, arnet::ArNet, msamples) where {T<:Integer}
    @extract arnet:H J p0 idxperm
    N = length(idxperm)
    length(x0) < N || error("Subsequence too long for the model")
    all(x -> 1 ≤ x ≤ 21, x0) || error("Subsequence numeric code should be in 1..21 ")
    l0 = length(x0)
    q = length(p0)
    N = length(H) # here N is N-1 !!
    backorder = sortperm(idxperm)

    W = Vector{Float64}(undef, msamples)
    res = Matrix{Int}(undef, N + 1, msamples)
    Threads.@threads for i in 1:msamples
        totH = Vector{Float64}(undef, q)
        sample_z = -ones(Int, N + 1)
        for k in 1:l0
            sample_z[backorder[k]] = x0[k]
        end

        if sample_z[1] == -1
            sample_z[1] = wsample(1:q, p0)
        end
        logw = log(p0[sample_z[1]])
        for site in 1:N
            Js = J[site]
            h = H[site]
            copy!(totH, h)
            @avx for i in 1:site
                for a in 1:q
                    totH[a] += Js[a, sample_z[i], i]
                end
            end
            p = softmax(totH)
            if sample_z[site+1] == -1
                sample_z[site+1] = wsample(1:q, p)
            end
            logw += log(p[sample_z[site+1]])
        end
        W[i] = exp(logw)
        res[:, i] .= sample_z
    end
    return W, permuterow!(res, backorder)
end

function permuterow!(x::AbstractMatrix, p::Vector)
    isperm(p) || error("not a permutation")
    for j in 1:size(x, 2)
        vx = @view x[:, j]
        Base.permute!(vx, p)
    end
    x
end

function permuterow!(x::AbstractVector, p::Vector)
    isperm(p) || error("not a permutation")
    Base.permute!(x, p)
end

"""
    loglikelihood(x0::Vector{T}, arnet::ArNet) where {T<:Integer})
Return the loglikelihood of sequence `x0` encoded in integer values in `1:q` under the model `arnet``. 
"""
(loglikelihood(x0::Vector{T}, arnet::ArNet) where {T<:Integer}) = sum(log, arnet(x0))

"""
    loglikelihood(x0::String, arnet::ArNet) where {T<:Integer})
Return the loglikelihood of the `String` `x0` under the model `arnet`. 
"""
function loglikelihood(s0::String, arnet::ArNet)
    x0 = letter2num.(collect(s0))
    sum(log, arnet(x0))
end

"""
    loglikelihood(x0::Matrix{T}, arnet::ArNet) where {T<:Integer}) 
Return the vector of loglikelihoods computed from `Matrix` `x0` under the model
`arnet`. `size(x0) == N,M` where `N` is the sequences length, and `M` the number
of sequences. The returned vector has `M` elements.
"""
(loglikelihood(x0::Matrix{T}, arnet::ArNet) where {T<:Integer}) = sum(log, arnet(x0), dims=1)[:]


# warning the gauge of H[:,1] is to be determined !!
function tensorize(arnet::ArNet; tiny::Float64=1e-16)
    @extract arnet:J H idxperm p0
    N = length(idxperm)
    q = length(H[1])
    p0pc = (1.0 - tiny) * p0 .+ tiny / q
    outJ = zeros(q, q, N, N)
    outH = zeros(q, N)
    shiftH0 = sum(log.(p0pc)) / q
    outH[:, idxperm[1]] .= log.(p0pc) .- shiftH0
    for i in 1:N-1
        si = idxperm[i+1]
        Js = J[i]
        outH[:, si] .= H[i]
        for j in 1:i
            sj = idxperm[j]
            outJ[:, :, si, sj] .= Js[:, :, j]
            outJ[:, :, sj, si] .= Js[:, :, j]'
        end
    end
    outJ, outH
end

"""
    dms_single_site(arnet::ArNet, arvar::ArVar, seqid::Int; pc::Float64=0.1)
    
Return a `q×L` matrix of containing `-log(P(mut))/log(P(seq))` for all single
site mutants of the reference sequence `seqid`, and a vector of the indices of
the residues of the reference sequence that contain gaps (i.e. the 21
amino-acid) for which the score has no sense and is set by convention to `+Inf`.
A negative value indicate a beneficial mutation, a value 0 indicate
the wild-type amino-acid.
"""
function dms_single_site(arnet::ArNet, arvar::ArVar, seqid::Int; pc::Float64=0.1)
    @extract arnet:H J p0 idxperm
    @extract arvar:Z M N q

    1 ≤ seqid ≤ M || error("seqid=$seqid should be in the interval [1,...,$M]")

    ppc = (1 - pc) * p0 + pc * ones(q) / q
    Da = fill(Inf64, q, N)
    xori = Z[:, seqid]
    xmut = copy(xori)
    idxnogap = findall(x -> x != 21, xori)
    arlike = zeros(N)
    arlike0 = zeros(N)
    _outputarnet!(arlike0, xori, J, H, ppc, N, q)
    ll0 = -sum(log.(arlike0))

    @inbounds for i in idxnogap
        if xori[i] == 21
            continue
        end
        for a in 1:q
            if a != xori[i]
                xmut[i] = a
                _outputarnet!(arlike, xmut, J, H, ppc, N, q)
                Da[a, i] = -sum(log.(arlike)) - ll0
            else
                Da[a, i] = 0.0
            end
        end
        xmut[i] = xori[i] #reset xmut to the original velue 
    end
    invperm = sortperm(idxperm)
    return Da[:, invperm], sort!(idxperm[setdiff(1:N, idxnogap)])
end
