module TestSchur

using Test
import JointDiag

include("utils.jl")

function test_lapack_exception()
    A = [
        -1.96262e-16 0 0 0
        -1.17757e-16 0 0 0
        0 0 0 0
        0 0 0 0
    ]
    B = [
        0 0 0 -3.92523e-17
        1 0 0 -1.17757e-16
        0 1 0 0
        0 0 1 0
    ]
    λ = [0.285173013664907, 0.714826986335093]
    # Throws `LAPACKException(1)`
    el = JointDiag._solve([A, B], λ, schur_solver)[2]
    testelements(el, reshape([0; 0], 2, 1); atol = 1e-10)
    return
end

function _test_cgt96_e51(solver)
    ɛ = 1e-4
    Iɛ = [
        1-ɛ 0
        0 1+ɛ
    ]
    J = [
        0 1
        1 0
    ]
    Z = zeros(2, 2)
    A = [
        Iɛ Z
        Z J
    ]
    B = [
        J Z
        Z Iɛ
    ]
    α = 0.219
    testelements(
        JointDiag._solve([A, B], [α, 1 - α], solver)[2],
        [1.0  -1.0; 1.0 1.0; -1.0 1.0]';
        rtol = 1e-7,
    )
    return
end

# Example 5.1 of CGT97
function test_cgt96_e51()
    @testset "Schur" begin
        _test_cgt96_e51(schur_solver)
    end
end

function runtests(args...)
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)(args...)
            end
        end
    end
end

end

TestSchur.runtests()
