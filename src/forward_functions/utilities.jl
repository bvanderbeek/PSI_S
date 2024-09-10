# -- utilities, do not touch >:(
function get_field(name,fields,rays2nodes,rays_outdom,nray,model)
    # for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray])
    #     nucleus = model.fields[model.fieldslist[name]].nodes2nuclei[j]
    #     fields[model.fieldslist[name],i] = model.fields[model.fieldslist[name]].v[nucleus]
    # end

    @. fields[model.fieldslist[name],1:rays2nodes[2,nray]-rays2nodes[1,nray]+1] = model.fields[model.fieldslist[name]].v[model.fields[model.fieldslist[name]].nodes2nuclei[rays2nodes[1,nray]:rays2nodes[2,nray]]]

    if rays_outdom[1 + 2*(model.fieldslist[name]-1),nray] != 0
        n1, n2 = rays_outdom[1 + 2*(model.fieldslist[name]-1),nray], rays_outdom[2 + 2*(model.fieldslist[name]-1),nray]
        fields[model.fieldslist[name],n1:n2] .= model.fields[model.fieldslist[name]].ref_value
    end
    return @view fields[model.fieldslist[name],begin:rays2nodes[2,nray]-rays2nodes[1,nray]+1]
end

function average_velocity(raytmp,rays2nodes,nray)
    np = rays2nodes[2,nray] - rays2nodes[1,nray] + 1
    @inbounds for i in 1:np
        if i != np
            raytmp.u_path[i] = raytmp.u[i]
        end
        if i != 1
            raytmp.u_path[i-1] += raytmp.u[i]
            raytmp.u_path[i-1] /= 2.0
        end
    end
end