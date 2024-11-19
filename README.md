# DecayTreeDataFrames

`DecayTreeDataFrames` is a Julia package for creating and manipulating tree structures, particularly useful for particle decay processes.

## Features

- **Tree Construction**: Build decay trees from nested tuples.
- **Data Propagation**: Propagate information up or down the tree structure.
- **Custom Transformations**: Apply transformations for physics-specific workflows.
- **Lorentz Transformations**: Handle and label transformations systematically.

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
transform_from_childen!(df, :subsystem) do logistics
    @unpack left, right = logistics
    i1, i2 = df[left, :subsystem], df[right, :subsystem]
	vcat(i1, i2)
end
```

###

Tuple notations for the node are also useful. They need to be collected from bottom to top, using `transform_from_childen!` function
```julia
df.notation = getproperty.(df.logistics, :info)
transform_from_childen!(df, :notation) do logistics
    @unpack left, right = logistics
    i1, i2 = df[left, :notation], df[right, :notation]
    "($i1, $i2)"
end
df
```

### Lorentz Transformations
Apply and propagate Lorentz transformations:
```julia
df.transformation .= "I"
transform_from_parent!(f, df, :transformation)
```

## License

This package is released under the MIT License.