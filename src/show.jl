
function DD.show_after(io::IO, mime::MIME"text/plain", A::AbstractGeoArray)
    if parent(A) isa AbstractDiskArray 
        if parent(A) isa FileArray 
            printstyled(io, "\nfrom file:\n"; color=:light_black)
            print(io, filename(parent(A)))
        end
    else
        DD.show_array(io, mime, parent(A))
    end
end

function DD.show_after(io, mime, stack::AbstractGeoStack) 
    if !(window(stack) isa Nothing)
        printstyled(io, "\nwith window:\n"; color=:light_black)
        for (i, dim) in enumerate(window(stack))
            print(io, ' ')
            show(IOContext(io, :compact=>true), mime, dim)
            i != length(window(stack)) && print(io, ',')
        end
    end
    if data(stack) isa FileStack 
        printstyled(io, "\nfrom file:\n"; color=:light_black)
        println(io, filename(stack))
    end
end

# Stack types can be enourmous. Just use nameof(T)
function Base.summary(io::IO, ser::AbstractGeoSeries{T,N}) where {T,N}
    if N == 0  
        print(io, "0-dimensional ")
    elseif N == 1
        print(io, size(ser, 1), "-element ")
    else
        print(io, join(size(ser), "×"), " ")
    end
    printstyled(io, string(nameof(typeof(ser)), "{$(nameof(T)),$N}"); color=:blue)
end

function Base.show(io::IO, mime::MIME"text/plain", ser::AbstractGeoSeries{T,N}) where {T,N}
    lines = _print_array_info(io, mime, ser)
    ds = displaysize(io)
    ioctx = IOContext(io, :displaysize => (ds[1] - lines, ds[2]))
end

function Base.show(io::IO, mime::MIME"text/plain", mode::AbstractProjected)
    DD._printmode(io, mode)
    print(io, " - ")
    DD._printorder(io, mode)
    print(io, " ")
    DD._printspan(io, mode)
    print(io, " ")
    DD._printsampling(io, mode)
    if !(crs(mode) isa Nothing)
        print(io, " crs: ", nameof(typeof(crs(mode))))
    end
    if !(mappedcrs(mode) isa Nothing)
        print(io, " mappedcrs: ", nameof(typeof(mappedcrs(mode))))
    end
end

function _print_array_info(io, mime, A)
    lines = 0
    summary(io, A)
    DD._printname(io, name(A))
    lines += DD._printdims(io, mime, dims(A))
    !(isempty(dims(A)) || isempty(refdims(A))) && println(io)
    lines += DD._printrefdims(io, mime, refdims(A))
    println(io)
    return lines
end
