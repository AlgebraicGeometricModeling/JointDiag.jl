all: JointDiag.jl
	julia --color=yes make.jl

JointDiag.jl:
	git clone -b gh-pages https://github.com/AlgebraicGeometricModeling/JointDiag.jl.git

publish:
	cd JointDiag.jl; git add docs; git commit -am "update gh-pages"; git push
