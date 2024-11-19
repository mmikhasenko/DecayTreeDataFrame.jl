"""
    transform_tree_down!(func, df, entry_index)

Modify the data frame starting from the root node following to the rows of children.
The function `func` is applied to the row structure of the data frame.
It uses `df.logistics::Routing` to navigate.
The iteration over the tree starts from the the `entry_index` node.
"""
function transform_tree_down!(func, df, entry_index)
    @unpack left, right, me = df.logistics[entry_index]
    #
    func(df[me, :])
    left != 0 && transform_tree_down!(func, df, left)
    right != 0 && transform_tree_down!(func, df, right)
    return
end

"""
    transform_tree_down!(func, df, entry_index)

Modify the data frame starting from children and reaching the root node.
The function `func` is applied to the row structure of the data frame.
It uses `df.logistics::Routing` to navigate.
The iteration over the tree starts from the the `entry_index` node.
"""
function transform_tree_up!(func, tree, entry_index)
    @unpack left, right, me = tree.logistics[entry_index]
    #
    left != 0 && transform_tree_up!(func, tree, left)
    right != 0 && transform_tree_up!(func, tree, right)
    func(tree[me, :])
    return
end

"""
    transform_from_parent!(func, df, name)

Modify the data frame starting from children and reaching the root node.
The function `func` is applied the `logistic::Routing` object. The function is not applied to the root note.
"""
function transform_from_parent!(func, df, name)
    transform_tree_down!(df, size(df, 1)) do row
        @unpack parent, me = row.logistics
        value = parent != 0 ? func(row.logistics) : df[me, name]
        setproperty!(row, name, value)
    end
end

"""
    transform_from_childen!(func, df, name)

Modify the data frame starting from children and reaching the root node.
The function `func` is applied the `logistic::Routing` object. The function is not applied to the leaves of the tree (no children).
"""
function transform_from_children!(func, df, name)
    transform_tree_up!(df, size(df, 1)) do row
        @unpack left, right, me = row.logistics
        value = complex(left, right) != 0 + 0im ? func(row.logistics) : df[me, name]
        setproperty!(row, name, value)
    end
    df
end
