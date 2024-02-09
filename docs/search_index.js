var documenterSearchIndex = {"docs":
[{"location":"expl/1.joint_diag/#An-example-with-simple-eigenvalues","page":"An example with simple eigenvalues","title":"An example with simple eigenvalues","text":"","category":"section"},{"location":"expl/1.joint_diag/","page":"An example with simple eigenvalues","title":"An example with simple eigenvalues","text":"using JointDiag, LinearAlgebra\n\nN = 4\nE = randn(N,N)\n\nXi = [1 1 1 1   ;\n\t  1 -1 2 2  ;\n\t  -1 -1 0 3 ]\n\t  \nM00 = E*diagm(Xi[1,:])*inv(E)\nM01 = E*diagm(Xi[2,:])*inv(E)\nM02 = E*diagm(Xi[3,:])*inv(E)\n\nM0 = [M00, M01, M02]\n\nX0, E0, Slv = joint_diag(M0, NewtonJointDiag())\n","category":"page"},{"location":"expl/1.joint_diag/","page":"An example with simple eigenvalues","title":"An example with simple eigenvalues","text":"The joint eigenvalue vectors are:","category":"page"},{"location":"expl/1.joint_diag/","page":"An example with simple eigenvalues","title":"An example with simple eigenvalues","text":"julia> X0\n3×4 Matrix{Float64}:\n 1.0   1.0           1.0   1.0\n 2.0   2.0          -1.0   1.0\n 3.0  -3.03672e-17  -1.0  -1.0","category":"page"},{"location":"expl/1.joint_diag/","page":"An example with simple eigenvalues","title":"An example with simple eigenvalues","text":"The common eigenvectors are:","category":"page"},{"location":"expl/1.joint_diag/","page":"An example with simple eigenvalues","title":"An example with simple eigenvalues","text":"julia> E0\n4×4 Matrix{Float64}:\n -0.42709    0.785318     0.18       -0.519269\n -0.824614  -0.45775     -0.0883551  -0.482747\n -0.229514   0.00303585   0.431565    0.588898\n -0.291424  -0.41681      0.879514    0.387962","category":"page"},{"location":"expl/1.joint_diag/","page":"An example with simple eigenvalues","title":"An example with simple eigenvalues","text":"We verify that the residual error of the decomposition is small:","category":"page"},{"location":"expl/1.joint_diag/","page":"An example with simple eigenvalues","title":"An example with simple eigenvalues","text":"julia> norm.(M0-[E0*diagm(X0[i,:])*inv(E0) for i in 1:length(M0)])\n3-element Vector{Float64}:\n 6.926513892003204e-16\n 1.0466060440537818e-14\n 3.224070268441243e-15","category":"page"},{"location":"code/joint_diag/#Joint-Diagonalization","page":"Joint Diagonalization","title":"Joint Diagonalization","text":"","category":"section"},{"location":"code/joint_diag/","page":"Joint Diagonalization","title":"Joint Diagonalization","text":"Pages = [\"joint_diag.md\"]","category":"page"},{"location":"code/joint_diag/","page":"Joint Diagonalization","title":"Joint Diagonalization","text":"JointDiag.joint_diag\nJointDiag.joint_reduce\nJointDiag.JRS\nJointDiag.ReorderedSchurSolver\nJointDiag.AbstractSolver\nJointDiag.cluster_eigenvalues","category":"page"},{"location":"code/joint_diag/#JointDiag.joint_diag","page":"Joint Diagonalization","title":"JointDiag.joint_diag","text":" joint_diag(M::Vector{Matrix{C}}, Solver::RandJointDiag)\n\nCompute the joint diagonalization of an array M of square matrices M[1],...,M[n] using a random combination of the matrices and its Schur factorization to get the common eigenvectors. It outputs\n\nX the vectors of eigenvalues, which are the columns of X.\nE the common eigenvectors such that M[i]*E=E*diagm(X[i,:])\n\n\n\n\n\n joint_diag(M::Vector{Matrix{C}}, Solver::NewtonJointDiag)\n\nCompute the joint diagonalization of an array M of square matrices M[1], ..., M[n] by applying a Newton-type iteration to minimize the off-diagonal norm of the matrices F*M[i]*E with the constraint F*E=I. It outputs\n\nX the vectors of eigenvalues, which are the columns of X.\nE the common eigenvectors such that M[i]*E=E*diagm(X[i,:])\n\nIt implements the method described in\n\n[KMY22] Khouja, Rima, Mourrain, Bernard, and Yakoubsohn, Jean-Claude. Newton-type methods for simultaneous matrix diagonalization. Calcolo 59.4 (2022): 38. doi:10.1007/s10092-022-00484-3, https://hal.science/hal-03390265.\n\n\n\n\n\n joint_diag(M::Vector{Matrix{C}}, Solver::JacobiJointDiag)\n\nCompute the joint diagonalization of an array M of square matrices M[1], ..., M[n] by applying a Jacobi-type iteration to minimize the off-diagonal norm of the matrices F*M[i]*E with the constraint F*E=I. It outputs\n\nX the vectors of eigenvalues, which are the columns of X.\nE the common eigenvectors such that M[i]*E=E*diagm(X[i,:])\n\n\n\n\n\n","category":"function"},{"location":"code/joint_diag/#JointDiag.joint_reduce","page":"Joint Diagonalization","title":"JointDiag.joint_reduce","text":" joint_reduce(H::Vector{Matrix{C}}, Solver = JRS())\n\nCompute the joint diagonal reduction of an array H of square matrices H[1],...,H[n] using a random combination of the matrices, SVD and rank reduction to obtain vectors of eigenvalues X,  left factor U, right factor V such that H=[U*diagm(X[i,:])*V for i in 1:length(H)] It outputs\n\nX the eigenvalue matrix, which columns are the vectors of common eigenvalues\nU the left factor, which columns are  orthogonal\nV the right factor, which rows are orthogonal\n\n\n\n\n\n","category":"function"},{"location":"code/joint_diag/#JointDiag.JRS","page":"Joint Diagonalization","title":"JointDiag.JRS","text":"Describe the Joint Reduce Solver\n\n - `rank` is a rank function used to determine the rank from the singular values\n - `diag_solver`is the solver used for the joint diagonalization\n\n\n\n\n\n","category":"type"},{"location":"code/joint_diag/#JointDiag.ReorderedSchurSolver","page":"Joint Diagonalization","title":"JointDiag.ReorderedSchurSolver","text":"struct ReorderedSchurSolver{T,RNGT<:Random.AbstractRNG} <: AbstractSolver\n    ɛ::T\n    rng::RNGT\nend\n\nSimultaneous diagonalization of commuting matrices using the method of [CGT97].\n\n[CGT97] Corless, R. M.; Gianni, P. M. & Trager, B. M. A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots Proceedings of the 1997 international symposium on Symbolic and algebraic computation, 1997, 133-140\n\n\n\n\n\n","category":"type"},{"location":"code/joint_diag/#JointDiag.AbstractSolver","page":"Joint Diagonalization","title":"JointDiag.AbstractSolver","text":"AbstractSolver\n\nSolver for joint diagonalization of matrices.\n\n\n\n\n\n","category":"type"},{"location":"code/joint_diag/#JointDiag.cluster_eigenvalues","page":"Joint Diagonalization","title":"JointDiag.cluster_eigenvalues","text":"cluster_eigenvalues(_atol, v)\n\nClustering the values v following [CGT97].\n\n[CGT97] Corless, R. M.; Gianni, P. M. & Trager, B. M. A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots Proceedings of the 1997 international symposium on Symbolic and algebraic computation, 1997, 133-140\n\n\n\n\n\n","category":"function"},{"location":"#JointDiag","page":"Home","title":"JointDiag","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Package for joint diagonalisation of pencils of matrices","category":"page"}]
}