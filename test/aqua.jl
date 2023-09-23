import Aqua
import IsingModels
using Test: @testset

@testset "aqua" begin
    Aqua.test_all(IsingModels; ambiguities = false)
end
