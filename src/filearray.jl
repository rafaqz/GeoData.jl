
struct FileArray{X,T,N,F<:AbstractString,S<:Tuple,K,EC,HC} <: AbstractDiskArray{T,N}
    filename::F
    size::S
    key::K
    eachchunk::EC
    haschunks::HC
    write::Bool
end
function FileArray{X,T,N}(
    filename::F, size::S, key::K, eachchunk::EC=nothing, 
    haschunks::HC=DA.Unchunked(), write=false
) where {X,T,N,F,S<:Tuple,K,EC,HC}
    FileArray{X,T,N,F,S,K,EC,HC}(filename, size, key, eachchunk, haschunks, write)
end
function FileArray{X,T,N}(filename::AbstractString, size::Tuple; 
    key=nothing, eachchunk=size, haschunks=DA.Unchunked(), write=false
) where {X,T,N}
    FileArray{X,T,N}(filename, size, key, eachchunk, haschunks, write)
end

ConstructionBase.constructorof(::Type{<:FileArray{X,T,N}}) where {X,T,N} = FileArray{X,T,N}

filename(A::FileArray) = A.filename
key(A::FileArray) = A.key
Base.size(A::FileArray) = A.size
DA.eachchunk(A::FileArray) = A.eachchunk
DA.haschunks(A::FileArray) = A.haschunks

function Base.open(f::Function, A::FileArray{X}; write=A.write, kw...) where X
    _read(f, X, filename(A); key=key(A), write, kw...)
end


DA.readblock!(A::FileArray, dst, r::AbstractUnitRange...) = open(o -> dst .= o[r...], A)
DA.writeblock!(A::FileArray, src, r::AbstractUnitRange...) = 
    open(o -> o[r...] .= src, A; write=A.write)


struct GeoDiskArray{T,N,V<:AbstractArray{T,N},EC,HC} <: AbstractDiskArray{T,N}
    var::V
    eachchunk::EC
    haschunks::HC
end
GeoDiskArray(var) = GeoDiskArray(var, _eachchunk(var), _haschunks(var))

Base.parent(A::GeoDiskArray) = A.var
Base.size(A::GeoDiskArray{T,N}) where {T,N} = size(parent(A))::NTuple{N,Int}

DA.haschunks(A::GeoDiskArray) = A.haschunks
DA.eachchunk(A::GeoDiskArray) = A.eachchunk
DA.readblock!(A::GeoDiskArray, aout, r::AbstractUnitRange...) = aout .= parent(A)[r...]
DA.writeblock!(A::GeoDiskArray, v, r::AbstractUnitRange...) = parent(A)[r...] .= v
