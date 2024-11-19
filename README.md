# DecayTreeDataFrames

`DecayTreeDataFrames` is a Julia package for creating and manipulating tree structures using a linear representation of a tree. A table interface is used for the storing the information. A navigation through the tree is done by using a `Rooting` structure, that just contains row indices of the table for `left` (child 1), `right` (child 2) and `parent`, and `me` (this) nodes.

## Installation

Include the module in your project:
```julia
using DecayTreeDataFrames
```

## Usage Example

### Build a Tree
```julia
tree_structure = ((:e, (:r, :u)), (:h, :g))
logistics = build_tree(tree_structure)
df = DataFrame(; logistics)
```

### Content of subsystem

For every node, one need to collect the full list of leaves. This can be done by:
```julia
df.subsystem = getproperty.(df.logistics, :info) .|> vcat
transform_from_children!(df, :subsystem) do logistics
    @unpack left, right = logistics
    i1, i2 = df[left, :subsystem], df[right, :subsystem]
    vcat(i1, i2)
end
```

### Tuple notations

Tuple notations for the node are also useful. They need to be collected from bottom to top, using `transform_from_children!` function
```julia
df.notation = getproperty.(df.logistics, :info)
transform_from_children!(df, :notation) do logistics
    @unpack left, right = logistics
    i1, i2 = df[left, :notation], df[right, :notation]
    "($i1, $i2)"
end
```

### Transition chains

For every node,
we can generate a string that shows how to arrive to this node from the root. This can be done by:
```julia
df.chain .= "0"
transform_from_parent!(df, :chain) do logistics
    @unpack parent, me, left, right = logistics
    df[parent, :chain] * " -> " * string(df[me, :logistics].info)
end
```

The call `transform_from_parent!` does not apply the function to the root node.
While `transform_from_children!` does not apply the function to the leaves.


## License

This package is released under the MIT License.