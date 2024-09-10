function delayed_rejection(model,update,observables,vnox,nodes2rays,rays2nodes,rays_outdom,IP,ind,action)
    update.trial_sample.n .= 2
    if action == 1
        pertv(model.fields[update.field[1]],IP.MCS.pert_size/IP.DR.rescale,vnox,nodes2rays,rays2nodes,update,IP;ind=ind)
    elseif action == 2
        pertp(model.fields[update.field[1]],IP.MCS.pert_size/IP.DR.rescale,vnox,nodes2rays,rays2nodes,update,IP;ind=ind)
    end
    update_predictions(model,observables,vnox,rays2nodes,rays_outdom,update)
    evaluate_model(model,observables,update,IP,vnox,nodes2rays,rays2nodes,rays_outdom)
end

function acceptance_ratio_delayed_rejection(model,update,IP)
    voronoi = model.fields[update.field[1]]
    p = IP.MCS.Lnorm
    q_ratio = 0.0
    if update.action[1] == 1
        dv1, dv2 = update.trial_sample.dv1, update.trial_sample.dv2
        σ = IP.MCS.pert_size * (voronoi.vlims[2]-voronoi.vlims[1])
        if (σ != 0)
            q_ratio = -0.5*((dv1[1]-dv2[1])^2/σ^2) + 0.5*(dv1[1]^2/σ^2)
        end
    elseif update.action[1] == 2
        ds1, ds2 = update.trial_sample.ds1, update.trial_sample.ds2
        q_ratio = 0.0
        for i in eachindex(ds1)
            σ = IP.MCS.pert_size * (voronoi.slims[i][2]-voronoi.slims[i][1])
            (σ == 0) && continue
            q_ratio += -0.5*(((ds1[i]-ds2[i])^2)/σ^2) + 0.5*(ds1[i]^2/σ^2)
        end
    end
    if update.term[1] == -Inf
        return -Inf
    else
        α1 = min(1,exp(-(1/p)*(update.trial_sample.misfit[1] - update.misfit[1])/model.T[1])*exp(update.trial_sample.term[1]))
        α2 = min(1,exp(-(1/p)*(update.trial_sample.misfit[1] - model.misfit[1])/model.T[1])*exp(update.trial_sample.term[1]))
        β1 = (1-α1)
        β2 = (1-α2)
        acceptance_ratio = -(1/p)*(update.misfit[1]-model.misfit[1])/model.T[1] + q_ratio + log(β1/β2)
        # print("\nacc",acceptance_ratio," ",update.misfit[1]-model.misfit[1]," ",q_ratio," ",log(β1/β2))
        return acceptance_ratio
    end
end