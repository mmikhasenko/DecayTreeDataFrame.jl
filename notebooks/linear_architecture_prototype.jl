### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 75b6b4d3-5a64-4f1e-a40d-1b476f8b6b7d
begin
    using Setfield
    using DataFrames
    using Parameters
end

# ╔═╡ a04ac952-4289-42a0-a693-e2251473b225
md"""
# Prototype for linear architecture

Idea was to store the tree nodes in an array where the children are registered with array index.
"""

# ╔═╡ dac292e3-a19a-4ddf-9c4b-6c08949ce524
struct Routing{T}
    me::Int
    left::Int  # particle one
    right::Int # particle two
    parent::Int
    info::T
end

# ╔═╡ b168da5c-10fd-417b-9393-8928e60cacea
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

# ╔═╡ 526a211c-456f-428b-b64d-97f49cbc9519
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

# ╔═╡ 471d2e67-868a-4fd3-98db-af890b96bbd6
function build_tree(tuple)
    _nodes = []
    root_index = build_tree(tuple, _nodes)
    # add parents
    _parents = parents(_nodes)
    map_indices = map(enumerate(zip(_parents, _nodes))) do (i, (ip, node))
        Routing(i, node.left, node.right, ip, node.info)
    end
end

# ╔═╡ 9432334c-d2c1-4b0a-a861-f1610f245b38
function modify_tree_down(f!, tree, entry_index)
    @unpack left, right, parent, me = tree.logistics[entry_index]
    #
    f!(tree[me, :])
    left != 0 && modify_tree_down(f!, tree, left)
    right != 0 && modify_tree_down(f!, tree, right)
    return
end

# ╔═╡ 7b74d569-fd81-4296-a8d5-dae029726a7a
function modify_tree_up(f!, tree, entry_index)
    @unpack left, right, parent, me = tree.logistics[entry_index]
    #
    left != 0 && modify_tree_up(f!, tree, left)
    right != 0 && modify_tree_up(f!, tree, right)
    f!(tree[me, :])
    return
end

# ╔═╡ ae211377-05b3-4d6a-b819-b42070c508d3
function propagate_from_parent!(f, df, name)
    modify_tree_down(df, nrow(df)) do row
        @unpack parent, me = row.logistics
        value = parent != 0 ? f(row.logistics) : df[me, name]
        setproperty!(row, name, value)
    end
end

# ╔═╡ 712b0390-27fb-4cd0-ab76-8d1bdcee9f09
function propagate_from_childen!(f, df, name)
    modify_tree_up(df, nrow(df)) do row
        @unpack left, right, me = row.logistics
        value = complex(left, right) != 0 + 0im ? f(row.logistics) : df[me, name]
        setproperty!(row, name, value)
    end
    df
end

# ╔═╡ b33f4073-434b-4f90-92ad-6d512cdc0913
md"""
## Application
"""

# ╔═╡ b4e8e8ea-b416-4876-84e3-35243487575b
tree_structure = ((:e, (:r, :u)), (:h, :g))

# ╔═╡ 162685ef-ec60-4bcf-aa93-7b1b6f6fa57c
logistics = build_tree(tree_structure)

# ╔═╡ 7d2e372d-2d98-4051-82a5-0a91e9611f04
df = DataFrame(; logistics)

# ╔═╡ a4a94f22-6b93-4241-b898-13ab1451236e
let
    df.notation = getproperty.(df.logistics, :info)
    propagate_from_childen!(df, :notation) do logistics
        @unpack left, right = logistics
        i1, i2 = df[left, :notation], df[right, :notation]
        "($i1, $i2)"
    end
    df
end

# ╔═╡ 3a630e49-45bd-4260-8bf4-a9db42d90ccf
let
    df.subsystem = getproperty.(df.logistics, :info) .|> vcat
    propagate_from_childen!(df, :subsystem) do logistics
        @unpack left, right = logistics
        i1, i2 = df[left, :subsystem], df[right, :subsystem]
        vcat(i1, i2)
    end
    df
end

# ╔═╡ 69810976-7abd-45cf-96f3-7a20c248d1b1
let
    df.chain .= "0"
    propagate_from_parent!(df, :chain) do logistics
        @unpack parent, me, left, right = logistics
        df[parent, :chain] * " -> " * string(df[me, :logistics].info)
    end
    df
end

# ╔═╡ 97e82d88-ba47-4a9a-8cc7-9a0a589f10eb
md"""
## Lorentz Transformations
"""

# ╔═╡ b73b65f7-e747-4c99-849f-3f5482f59e21
let
    df.transformation .= "I"
    propagate_from_parent!(df, :transformation) do logistics
        @unpack parent, me, left, right = logistics
        # 
        l_me = df[me, :logistics].info
        l_pe = df[parent, :logistics].info
        label = l_me * " -> " * l_pe
        # determine transformation: `particle-one` or `particle-two`
        which = df[parent, :logistics].left == me ? "Tf" : "Ts"
        my_transform = which * "⁻¹($label)"
        my_transform * " ∘ " * df[parent, :transformation]
    end
    df
end

# ╔═╡ 870aadb7-2054-493b-8fe1-2c98a454e365
let
    transform!(df, :transformation => ByRow() do t
        t * " [p1,p2, ...]"
    end => :four_vectors)
end

# ╔═╡ 432f56a5-475b-43a5-b3d4-404362ce85d9
let
    df.angles .= ""
    propagate_from_childen!(df, :angles) do logistics
        @unpack left = logistics
        s = df[left, :subsystem]
        "angles of p̄_$s in [p̄1,p̄2, ...]"
    end
    df
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
Setfield = "efcf1570-3423-57d1-acb7-fd33fddbac46"

[compat]
DataFrames = "~1.7.0"
Parameters = "~0.12.3"
Setfield = "~1.1.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "5fa5ec67d7902021cbdf3dc72d4b5a7a4e266c24"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "d0553ce4031a081cc42387a9b9c8441b7d99f32d"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.7"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╟─a04ac952-4289-42a0-a693-e2251473b225
# ╠═75b6b4d3-5a64-4f1e-a40d-1b476f8b6b7d
# ╠═dac292e3-a19a-4ddf-9c4b-6c08949ce524
# ╠═b168da5c-10fd-417b-9393-8928e60cacea
# ╠═526a211c-456f-428b-b64d-97f49cbc9519
# ╠═471d2e67-868a-4fd3-98db-af890b96bbd6
# ╠═9432334c-d2c1-4b0a-a861-f1610f245b38
# ╠═7b74d569-fd81-4296-a8d5-dae029726a7a
# ╠═ae211377-05b3-4d6a-b819-b42070c508d3
# ╠═712b0390-27fb-4cd0-ab76-8d1bdcee9f09
# ╟─b33f4073-434b-4f90-92ad-6d512cdc0913
# ╠═b4e8e8ea-b416-4876-84e3-35243487575b
# ╠═162685ef-ec60-4bcf-aa93-7b1b6f6fa57c
# ╠═7d2e372d-2d98-4051-82a5-0a91e9611f04
# ╠═a4a94f22-6b93-4241-b898-13ab1451236e
# ╠═3a630e49-45bd-4260-8bf4-a9db42d90ccf
# ╠═69810976-7abd-45cf-96f3-7a20c248d1b1
# ╟─97e82d88-ba47-4a9a-8cc7-9a0a589f10eb
# ╠═b73b65f7-e747-4c99-849f-3f5482f59e21
# ╠═870aadb7-2054-493b-8fe1-2c98a454e365
# ╠═432f56a5-475b-43a5-b3d4-404362ce85d9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
