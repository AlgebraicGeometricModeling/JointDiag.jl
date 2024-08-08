module TestSchur

using Test
import JointDiag

include("utils.jl")

function test_lapack_exception(T)
    A = T[
        -1.96262e-16 0 0 0
        -1.17757e-16 0 0 0
        0 0 0 0
        0 0 0 0
    ]
    B = T[
        0 0 0 -3.92523e-17
        1 0 0 -1.17757e-16
        0 1 0 0
        0 0 1 0
    ]
    λ = T[0.285173013664907, 0.714826986335093]
    # Throws `LAPACKException(1)`
    el = JointDiag.joint_diag([A, B], λ, schur_solver)[2]
    testelements(el, reshape(T[0; 0], 2, 1); atol = 1e-10)
    return
end

function _test_cgt96_e51(T, solver)
    ɛ = T(1e-4)
    Iɛ = T[
        1-ɛ 0
        0 1+ɛ
    ]
    J = T[
        0 1
        1 0
    ]
    Z = zeros(T, 2, 2)
    A = T[
        Iɛ Z
        Z J
    ]
    B = T[
        J Z
        Z Iɛ
    ]
    α = T(0.219)
    testelements(
        JointDiag.joint_diag([A, B], [α, 1 - α], solver)[2],
        T[1.0 -1.0; 1.0 1.0; -1.0 1.0]';
        rtol = 1e-7,
    )
    return
end

# Example 5.1 of CGT97
function test_cgt96_e51(T)
    @testset "Schur" begin
        _test_cgt96_e51(T, schur_solver)
    end
end

function test_empty(T)
    for n in 1:2
        testelements(
            JointDiag.joint_diag([zeros(T, 0, 0) for _ in 1:n], schur_solver)[2],
            zeros(T, n, 0),
        )
    end
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name) $T" for T in [Float64, BigFloat]
                getfield(@__MODULE__, name)(T)
            end
        end
    end
end

end

TestSchur.runtests()
