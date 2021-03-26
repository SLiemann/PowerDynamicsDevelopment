using ModelingToolkit

@variables t
@variables x[1:100](t)


@time
function myfun(x::Array{Num,1})
    a = Vector{Differential}()
    for i in x
        push!(a,Differential(i))
    end
    return a::Vector{Differential}
end
