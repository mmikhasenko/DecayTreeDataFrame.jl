using DecayTreeDataFrames
using Test
using DataFrames
using Parameters

tree_structure = ((:e, (:r, :u)), (:h, :g))
logistics = build_tree(tree_structure)

@test length(logistics) == 9
@test all(getproperty.(logistics, :me) .== 1:9)

df = DataFrame(; logistics)


df.notation = getproperty.(df.logistics, :info)
transform_from_children!(df, :notation) do logistics
    @unpack left, right = logistics
    i1, i2 = df[left, :notation], df[right, :notation]
    "($i1, $i2)"
end

@testset "Add notations" begin
    @test df.notation == [
        "e",
        "r",
        "u",
        "(r, u)",
        "(e, (r, u))",
        "h",
        "g",
        "(h, g)",
        "((e, (r, u)), (h, g))"]
end

df.subsystem = getproperty.(df.logistics, :info) .|> vcat
transform_from_children!(df, :subsystem) do logistics
    @unpack left, right = logistics
    i1, i2 = df[left, :subsystem], df[right, :subsystem]
    vcat(i1, i2)
end

@testset "Add subsystem content" begin
    @test df.subsystem == [
        ["e"],
        ["r"],
        ["u"],
        ["r", "u"],
        ["e", "r", "u"],
        ["h"],
        ["g"],
        ["h", "g"],
        ["e", "r", "u", "h", "g"]]
end


df.chain .= "0"
transform_from_parent!(df, :chain) do logistics
    @unpack parent, me, left, right = logistics
    df[parent, :chain] * " -> " * string(df[me, :logistics].info)
end

@testset "Add transition chains" begin
    @test df.chain == [
        "0 -> ⊡ -> e",
        "0 -> ⊡ -> ⊡ -> r",
        "0 -> ⊡ -> ⊡ -> u",
        "0 -> ⊡ -> ⊡",
        "0 -> ⊡",
        "0 -> ⊡ -> h",
        "0 -> ⊡ -> g",
        "0 -> ⊡",
        "0"]
end
