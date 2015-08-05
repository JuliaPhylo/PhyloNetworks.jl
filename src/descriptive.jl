# functions to describe a HybridNetwork object to avoid accessing the attributes directly
# Claudia August 2015

function tips(net::HybridNetwork)
    return [l.name for l in net.leaf]
end
