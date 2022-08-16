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
