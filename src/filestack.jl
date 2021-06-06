abstract type AbstractFileStack end

struct FileStack{X,K,F<:AbstractString,T<:NamedTuple,S<:NamedTuple,EC,HC} <: AbstractFileStack
    filename::F
    types::T
    sizes::S
    eachchunk::EC
    haschunks::HC
    write::Bool
end
function FileStack{X,K}(
    filename::F, types::T, sizes::S, eachchunk::EC, haschunks::HC, write::Bool
) where {X,K,F,T,S,EC,HC}
    FileStack{X,K,F,T,S,EC,HC}(filename, types, sizes, eachchunk, haschunks, write)
end

filename(fs::FileStack) = fs.filename
Base.keys(fs::FileStack{<:Any,K}) where K = K
Base.values(fs::FileStack{<:Any}) = (fs[k] for k in keys(fs))
function Base.getindex(fs::FileStack{X}, key::Symbol) where X
    size = fs.sizes[key]
    eachchunk = fs.eachchunk[key]
    haschunks = fs.haschunks[key]
    T = fs.types[key]
    N = length(size)
    FileArray{X,T,N}(filename(fs), size, key, eachchunk, haschunks, fs.write)
end
