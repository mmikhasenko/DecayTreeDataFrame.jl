struct Routing{T}
    me::Int
    left::Int  # particle one
    right::Int # particle two
    parent::Int
    info::T
end

function build_tree(tuple, nodes)
    if isa(tuple, Tuple)
        first_child = build_tree(tuple[1], nodes)
        second_child = build_tree(tuple[2], nodes)
        node = Routing(0, first_child, second_child, -1, "⊡")
        push!(nodes, node)
        return length(nodes)
    else
        node = Routing(0, 0, 0, -1, tuple |> string)
        push!(nodes, node)
        return length(nodes)
    end
end

function parents(tree)
    _pairs = map(enumerate(tree)) do (i, node)
        node.left => i, node.right => i
    end
    _list = vcat(getindex.(_pairs, 1), getindex.(_pairs, 2))
    only_nonzero = filter(x -> x[1] != 0, _list)
    sorted = sort(only_nonzero, by = x -> x[1])
    _parents = getindex.(sorted, 2)
    return append!(_parents, 0)
end

function build_tree(tuple)
    _nodes = []
    build_tree(tuple, _nodes)
    _parents = parents(_nodes)
    nodes = map(enumerate(zip(_parents, _nodes))) do (i, (ip, node))
        Routing(i, node.left, node.right, ip, node.info)
    end
    return nodes
end