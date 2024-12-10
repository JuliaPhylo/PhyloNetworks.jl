@testset "deprecated" begin
@test_throws "snaq! has been moved to the package SNaQ.jl" snaq!()
end