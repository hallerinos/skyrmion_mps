struct Bond
    s1::Int
    s2::Int
    r1::Vector{Float64}
    r2::Vector{Float64}
    desc::String
end

function Bond(s1::Int, s2::Int)
    return Bond(s1, s2, [0.0], [0.0], "")
end

Bond(s1::Int, s2::Int, r1::Vector, r2::Vector, desc::String="") = Bond(s1, s2, convert(Vector{Float64}, r1), convert(Vector{Float64}, r2), desc);

const Graph = Vector{Bond}