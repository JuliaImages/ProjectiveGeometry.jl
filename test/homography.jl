using Base.Test

@testset "homography" begin
    src=CartesianIndex{2}[]
    push!(src, CartesianIndex(1,1))
    push!(src, CartesianIndex(2,1))
    push!(src, CartesianIndex(3,1))
    push!(src, CartesianIndex(3,3))
    push!(src, CartesianIndex(1,3))

    des=CartesianIndex{2}[]
    push!(des, CartesianIndex(1,1))
    push!(des, CartesianIndex(3,1))
    push!(des, CartesianIndex(5,1))
    push!(des, CartesianIndex(4,4))
    push!(des, CartesianIndex(2,4))

    H = gethomography(src, des)

    @test H[3,3]==1

    p=H*[1, 1, 1]
    p/=p[3]
    @test_approx_eq p[1] 1
    @test_approx_eq p[2] 1

    p=H*[1, 2, 1]
    p/=p[3]
    @test_approx_eq p[1] 5/3
    @test_approx_eq p[2] 3
end