filter_ext(path, ext::AbstractString) = filter(fn -> splitext(fn)[2] == ext, readdir(path))
filter_ext(path, exts::Union{Tuple,AbstractArray}) = 
    filter(fn -> splitext(fn)[2] in exts, readdir(path))
filter_ext(path, ext::Nothing) = readdir(path)

# Check that array order matches expectation
checkarrayorder(A, order::Order) = map(d -> checkarrayorder(d, order), dims(A))
checkarrayorder(A, order::Tuple) = map(checkarrayorder, dims(A), order)
checkarrayorder(dim::Dimension, order::Order) =
    arrayorder(dim) == order || @warn "Array order for `$(DD.basetypeof(order))` is `$(arrayorder(dim))`, usually `$order`"

checkindexorder(A, order::Order) = map(d -> checkindexorder(d, order), dims(A))
checkindexorder(A, order::Tuple) = map(checkindexorder, dims(A), order)
checkindexorder(dim::Dimension, order::Order) =
    indexorder(dim) == order || @warn "Array order for `$(DD.basetypeof(order))` is `$(indexorder(dim))`, usually `$order`"

cleankeys(keys) = Tuple(map(Symbol, keys))
