## SETUP ##
using JLD
using PSI_S

o_dir = string(@__DIR__, "/output/", ARGS[1], "/")
burn_in = parse(Int64,ARGS[3])
chains = Vector{ModelConst}()
for k = 1:parse(Int64,ARGS[2])
   newchain = load(string(o_dir,ARGS[1],".",(k),".jld"),"MarkovChain")
   [push!(chains,model) for model in newchain[burn_in:end]]
end
save(string(o_dir,ARGS[1],"_inv.jld"),"MarkovChain",chains)
