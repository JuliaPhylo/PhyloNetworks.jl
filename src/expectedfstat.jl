function expectedf2matrix!(
    net::HybridNetwork;
    keepinternal::Bool=false,
    checkpreorder::Bool=true,
)
    m = descendenceweight(net; checkpreorder=checkpreorder)
    P = m.V # also m[:tips]
    #= idea: f2[i,j] = Ω[i,i] + Ω[j,j] -2 Ω[i,j] where
    Ω = m[:tips] * Diagonal(edgelengths) * transpose(m[:tips])
    is a linear function of edge lengths.
    express: f2[i,j] = (nn × nn × ne) edgeweight array * edgelength vector
    =#
    nn,ne = size(P) # number of nodes, number of edges
    edgeweight = zeros(eltype(P), nn,nn,ne)
    # todo: calculate these weights
    f2mat = zeros(eltype(P), nn,nn)
    for e in net.edge
        f2mat .+= e.length .* edgeweight[:,:,e.number]
    end
    return f2mat
end

# similarly: f3[x;i,j] = Ω[x,x] + Ω[i,j] - Ω[x,i] -2 Ω[x,j]
# and f4[i1,i2; i3,i4] =
