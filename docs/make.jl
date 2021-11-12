using Documenter, SquareIsingModel
import SquareIsingModel as Ising

makedocs(
    modules=[SquareIsingModel],
    sitename="SquareIsingModel.jl"
)

deploydocs(repo = "github.com/cossio/SquareIsingModel.jl.git")
