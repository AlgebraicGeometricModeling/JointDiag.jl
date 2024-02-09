using Documenter
using JointDiag

dir = "mrkd"
Expl = map(file -> joinpath("expl", file), filter(x ->endswith(x, "md"), readdir(dir*"/expl")))
Code = map(file -> joinpath("code", file), filter(x ->endswith(x, "md"), readdir(dir*"/code")))

makedocs(
         sitename = "JointDiag",
         authors = "B. Legat, B. Mourrain",
         modules = [JointDiag],
         build = "JointDiag.jl/docs",
         source = dir,
         pages = Any[
                     "Home" => "index.md",
                     "Examples" => Expl,
                     "Functions & types" => Code
                     ],
         repo =  Remotes.GitHub("AlgebraicGeometricModeling", "JointDiag.jl"),
         doctest = false
         )

#deploydocs(
#           repo = Remotes.GitHub("AlgebraicGeometricModeling", "JointDiag.jl"),
#           target = "site"
#           )

