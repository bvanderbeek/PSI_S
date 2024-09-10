function parallel_tempering(misfits,obs_lengths,obs_noise,temperatures,temp_ind,nchains,TMatrix,DMatrix)
    swaps = randperm(nchains)
    for i in 1:floor(Int64,nchains/2)
        ind1, ind2 = swaps[2*i-1], swaps[2*i]
        noise1 = @view obs_noise[ind1,:]
        noise2 = @view obs_noise[ind2,:]
        r = tempering_swap(
            misfits[ind1],misfits[ind2], noise1, noise2, obs_lengths, temperatures[ind1],temperatures[ind2]
            )
        TMatrix[temp_ind[ind1],temp_ind[ind2]] += r
        TMatrix[temp_ind[ind2],temp_ind[ind1]] += r
        DMatrix[temp_ind[ind1],temp_ind[ind2]] += 1
        DMatrix[temp_ind[ind2],temp_ind[ind1]] += 1
        if r == 1
            T1, T2 = temperatures[ind1], temperatures[ind2]
            temperatures[ind1], temperatures[ind2] = T2, T1
            printstyled(string("\nSWAP -> "); color = :blue)
            printstyled(string(ind1," - ",ind2,"; ",round(T1,digits=2)," - ",round(T2,digits=2),"; ",round(misfits[ind1],digits=2)," - ",round(misfits[ind2],digits=2),"\n"); color = :red)
            i1, i2 = temp_ind[ind1], temp_ind[ind2]
            temp_ind[ind1], temp_ind[ind2] = i2, i1
        end
    end
end

function tempering_swap(misfit1,misfit2,noise1,noise2,obs_lengths,T1,T2)
    r = log(rand())
    term = 0.0
    for i in eachindex(obs_lengths)
        term += obs_lengths[i]*(log(noise1[1,i])-log(noise2[1,i]))
    end
    if r < (1/T1 - 1/T2)*(0.5*(misfit1 - misfit2) + term)
        return 1.0
    end
    return 0.0
end
