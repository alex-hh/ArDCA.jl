"""
Compute contacts by computing deltadeltaE i.e. energy difference
between pair of mutations and single mutations

this is a single background approximation to the highly intractable

p(xi, xj)/p(xi)p(xj)

(note intractability is due to impossibility of computing any of the marginals)

It seems fairly surprising that this works from a single background.
The point is that we compute energy difference of all possible mutations,
so we are doing something very similar to what's done to infer contacts from
experimental fitness measurements.

Probably the pairwise parameterisation somehow ensures that the background isn't
too important.

pc is a pseudocount
"""
function epistatic_score(
    arnet::ArNet, seqnumeric::Vector{T}; pc::Float64=0.1, min_separation::Int=1
) where T <: Integer
    @extract arnet : H J p0 idxperm
    q = length(p0)
    N = length(idxperm)
    println(string("Alphabet size ", q, " num columns ", N))

    Da = zeros(q,N)
    Dab = zeros(q,q,N,N)

    xmut = [copy(seqnumeric) for _ in 1:Threads.nthreads()]
    arlike = [zeros(N) for _ in 1:Threads.nthreads()]
    ppc = (1-pc) * p0 + pc * ones(q)/q
    #_outputarnet!(arlike,xmut, J, H, ppc, N, q)
    # E0 =  sum(log.(arlike))
    # E0=0.0
    @inbounds for i in 1:N
        Threads.@threads for a in 1:q
            xmut[Threads.threadid()][i] = a
            _outputarnet!(arlike[Threads.threadid()],xmut[Threads.threadid()], J, H, ppc, N, q)
            Da[a,i] = -sum(log.(arlike[Threads.threadid()]))
            xmut[Threads.threadid()][i] = seqnumeric[i]
        end
    end  
    @inbounds for i in 1:N-1
        Threads.@threads for a in 1:q
            xmut[Threads.threadid()][i] = a
            for j in i+1:N
                 for b in 1:q
                    xmut[Threads.threadid()][i] = a
                    xmut[Threads.threadid()][j] = b
                    _outputarnet!(arlike[Threads.threadid()],xmut[Threads.threadid()],J,H,ppc,N,q)        
                    Dab[b,a,j,i] = -sum(log.(arlike[Threads.threadid()]))
                    xmut[Threads.threadid()][i] = seqnumeric[i]
                    xmut[Threads.threadid()][j] = seqnumeric[j]
                end
            end
        end
    end
    
    Jret = zeros(q,q,N,N)
    @inbounds for i in 1:N-1 
        for j in i+1:N
            for a in 1:q
                for b in 1:q 
                    Jret[b,a,j,i] = Dab[b,a,j,i] - Da[b,j] - Da[a,i]
                    #Jret[a,b,i,j] = Jret[b,a,j,i]
                end
            end
        end
    end
    Jzsg = zsg(Jret)
    FN = compute_APC(Jzsg, N, q)
    score = compute_ranking(FN, min_separation)
    
    permtuple=Tuple{Int,Int,Float64}[]
    sizehint!(permtuple,length(permtuple))
    for s in score
        si,sj,val= idxperm[s[1]],idxperm[s[2]],s[3]
        if si > sj 
            si,sj = sj,si
        end
        push!(permtuple,(si,sj,val))
    end
    return permtuple
end

function epistatic_score(arnet::ArNet, sequence::String; pc::Float64=0.1, min_separation::Int=1)
    @extract arnet: p0 idxperm
    # convert sequence to numeric (c.f. loglikelihood code)
    seqnumeric = letter2num.(collect(sequence))
    return epistatic_score(arnet, seqnumeric, pc=pc, min_separation=min_separation)
end

# https://stackoverflow.com/questions/38166522/julia-function-with-different-parameter-subtypes
function epistatic_score(
    arnet::ArNet, Z::Array{T,2}, seqid::Int; pc::Float64=0.1,min_separation::Int=1
) where T <: Integer
    # run on sequence at position seqid in alignment stored in arvar
    M = size(Z,1)
    1 ≤ seqid ≤ M || error("seqid=$seqid should be in the interval [1,...,$M]")
    xori = Z[:,seqid]
    return epistatic_score(arnet, xori, pc=pc, min_separation=min_separation)
end

function zsg(J::Array{Float64,4})
    q,q,N,N = size(J)
    Jzsg = zeros(q,q,N,N)
    @inbounds for i in 1:N-1
        for j in i+1:N
            Jzsg[:,:,j,i] .= J[:,:,j,i] - repeat(sum(J[:,:,j,i], dims=1)/q, q, 1) - repeat(sum(J[:,:,j,i], dims=2)/q, 1, q) .+ sum(J[:,:,j,i])/q^2
            #Jzsg[:,:,i,j] .= Jzsg[:,:,j,i]'
        end 
    end
    Jzsg
end


# If `Z` has size \$N × M\$ (i.e. \$M\$ sequences of length \$N\$), the resulting vector \$P_i\$ has
# length \$N (q-1)\$ and contains \$N\$ blocks (one for each residue position), each block containing
# the frequencies of the amino-acids, weighted according to `W`. The frequency of the last symbol,
# which usually represents the gap, is omitted and can be recovered by normalization. The resulting
# matrix \$P_{ij}\$ has size \$N (q-1) × N (q-1)\$ and it also has a block structure, with \$N × N\$
# blocks, one for each pair of residues (the last row and column of each block are omitted and can be
# recovered by normalization).

function invCov(
    Z::Array{Ti,2};
    pseudocount::Real = 0.8,
    θ = :auto,
    max_gap_fraction::Real = 0.9,
    remove_dups::Bool = false
    ) where Ti <: Integer

    if remove_dups
        Z, _ = remove_duplicate_sequences(Z)
    end
    N, M = size(Z)
    q = Int(maximum(Z))
    q ≥ 32 && error("parameter q=$q is too big (max 31 is allowed)")

    Pi_true, Pij_true, Meff, _ = compute_weighted_frequencies(Z, q, θ)

    Pi, Pij = add_pseudocount(Pi_true, Pij_true, Float64(pseudocount), q)

    C = compute_C(Pi, Pij)

    mJ = inv(cholesky(C))
    return mJ
end


function gDCA(
    filename::AbstractString;
    pseudocount::Real = 0.8,
    θ = :auto,
    max_gap_fraction::Real = 0.9,
    score::Symbol = :frob,
    min_separation::Integer = 5,
    remove_dups::Bool = false
)
    check_arguments(filename, pseudocount, θ, max_gap_fraction, score, min_separation)

    Z = read_fasta_alignment(filename, max_gap_fraction)
    q = Int(maximum(Z))

    mJ = invCov(Z, pseudocount=pseudocount, θ=θ, remove_dups=remove_dups)

    if score == :DI
        S = compute_DI_gauss(mJ, C, q)
    else
        S = compute_FN(mJ, q)
    end

    S = correct_APC(S)

    R = compute_ranking(S, min_separation)

    return mJ, R
end


function check_arguments(filename, pseudocount, θ, max_gap_fraction, score, min_separation)
    aerror(s) = throw(ArgumentError(s))
    0 <= pseudocount <= 1 ||
        aerror("invalid pseudocount value: $pseudocount (must be between 0 and 1)")
    θ == :auto || (θ isa Real && 0 <= θ <= 1) ||
        aerror("invalid θ value: $θ (must be either :auto, or a number between 0 and 1)")
    0 <= max_gap_fraction <= 1 ||
        aerror("invalid max_gap_fraction value: $max_gap_fraction (must be between 0 and 1)")
    score in [:DI, :frob] ||
        aerror("invalid score value: $score (must be either :DI or :frob)")
    min_separation >= 1 ||
        aerror("invalid min_separation value: $min_separation (must be >= 1)")
    isfile(filename) ||
        aerror("cannot open file $filename")

    return true
end


function printrank(io::IO, R::Vector{Tuple{Int,Int,Float64}})
    for I in R
        @printf(io, "%i %i %e\n", I[1], I[2], I[3])
    end
end
printrank(R::Vector{Tuple{Int,Int,Float64}}) = printrank(STDOUT, R)

printrank(outfile::AbstractString, R::Vector{Tuple{Int,Int,Float64}}) = open(f->printrank(f, R), outfile, "w")

compute_C(Pi::Vector{Float64}, Pij::Matrix{Float64}) = Pij - Pi * Pi'

function compute_APC(J::Array{Float64,4},N,q)
    FN = fill(0.0, N,N)
    @inbounds for i=1:N-1
        for j=i+1:N
            FN[j,i] = norm(J[1:q-1,1:q-1,j,i],2)
            FN[i,j] = FN[j,i]
        end
    end
    FN = correct_APC(FN)
    return FN
end

function correct_APC(S::Matrix)
    N = size(S, 1)
    Si = sum(S, dims=1)
    Sj = sum(S, dims=2)
    Sa = sum(S) * (1 - 1/N)

    S -= (Sj * Si) / Sa
    return S
end

function compute_ranking(S::Matrix{Float64}, min_separation::Int = 5)

    N = size(S, 1)
    R = Array{Tuple{Int,Int,Float64}}(undef, div((N-min_separation)*(N-min_separation+1), 2))
    counter = 0
    for i = 1:N-min_separation, j = i+min_separation:N
        counter += 1
        R[counter] = (i, j, S[j,i])
    end

    sort!(R, by=x->x[3], rev=true)
    return R

end
