using Test

include("schur.jl")
for file in readdir(joinpath(@__DIR__))
    if startswith(file, "test") && endswith(file, ".jl")
        @testset "$file" begin
            include(file)
        end
    end
end
