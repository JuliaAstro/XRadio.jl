"""
How a store's arrays are cut into chunks.

A chunk is the unit Zarr compresses and reads, so its shape decides what a
partial read costs. The default suits reading a Measurement Set a scan at a
time: every baseline, channel and polarization of a run of time steps.
"""

"""
    DEFAULT_CHUNK_BYTES

The uncompressed size, 16 MiB, that [`write`](@ref) fits the chunks of each
dataset to along the dimensions not given explicitly.
"""
const DEFAULT_CHUNK_BYTES = 16 * 2^20

# The order open dimensions are cut in. Time goes first, so a chunk is a run of
# whole time steps. Polarization is not cut unless named: it has at most four
# entries, varies fastest on disk, and cross-hand work reads it whole.
const _CUT_FIRST = (:time, :frequency, :baseline_id)
const _NEVER_CUT = (:polarization,)

struct ChunkRule
    fixed::Dict{Symbol, Int}
    bytes::Int
end

"""
    ChunkRule(chunks::NamedTuple, chunk_bytes, trees)

The chunking `write` applies, checked against the datasets `trees` hold: every
name in `chunks` must be a dimension of one of them, and every size positive.
"""
function ChunkRule(chunks::NamedTuple, chunk_bytes::Integer, trees)
    chunk_bytes > 0 || throw(ArgumentError("chunk_bytes must be positive; got $chunk_bytes"))
    for (name, n) in pairs(chunks)
        n isa Integer && n > 0 || throw(
            ArgumentError("the chunk size for `$name` must be a positive integer; got $(repr(n))")
        )
    end
    have = Set{Symbol}()
    for tree in trees
        _collect_dimnames!(have, tree)
    end
    unknown = setdiff(keys(chunks), have)
    isempty(unknown) || throw(
        ArgumentError(
            "`chunks` names " * join(repr.(collect(unknown)), ", ") * ", which " *
                (length(unknown) == 1 ? "is not a dimension" : "are not dimensions") *
                " of the data; its dimensions are " * join(repr.(sort!(collect(have))), ", ")
        )
    )
    return ChunkRule(Dict{Symbol, Int}(pairs(chunks)), Int(chunk_bytes))
end

function _collect_dimnames!(have, tree)
    for d in dims(tree)
        push!(have, msv4name(d))
    end
    foreach(sub -> _collect_dimnames!(have, sub), values(branches(tree)))
    return have
end

"""
    chunk_sizes(rule::ChunkRule, tree) -> Dict{Symbol, Int}

The chunk size along each dimension of the dataset `tree` (not its
sub-datasets): those `rule` fixes, and the rest cut from the dataset's largest
layer until one of its chunks fits `rule.bytes`. Every layer of the dataset
takes the same sizes, so the chunks of its visibilities, flags and weights
cover the same samples.
"""
function chunk_sizes(rule::ChunkRule, tree)
    sizes = Dict{Symbol, Int}()
    for d in dims(tree)
        name = msv4name(d)
        sizes[name] = min(get(rule.fixed, name, length(d)), length(d))
    end
    isempty(keys(tree)) && return sizes
    largest = argmax(k -> length(tree[k]) * _itemsize(tree[k]), collect(keys(tree)))
    array = tree[largest]
    names = collect(map(msv4name, dims(array)))
    open = [n for n in names if !haskey(rule.fixed, n) && n ∉ _NEVER_CUT]
    order = vcat(
        [n for n in _CUT_FIRST if n in open],
        [n for n in reverse(names) if n in open && n ∉ _CUT_FIRST],
    )
    for name in order
        bytes = _itemsize(array) * prod(n -> sizes[n], names; init = 1)
        bytes <= rule.bytes && break
        rest = bytes ÷ sizes[name]
        sizes[name] = max(1, rule.bytes ÷ rest)
    end
    return sizes
end

# Strings have no fixed size; they label other data and never decide its chunks.
_itemsize(array) = isbitstype(eltype(array)) ? sizeof(eltype(array)) : 0
