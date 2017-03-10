"""
gethomography estimates the 3x3 homography matrix from src to des points
"""

function gethomography(src::Array{CartesianIndex{2},1}, des::Array{CartesianIndex{2},1})
    if length(src)!=length(des)
        error("Different number of src and des points!")
    end

    n_points=length(src)

    if n_points<4
        error("Not enough points to estimate projection matrix")
    end

    xs = zeros(n_points)
    ys = zeros(n_points)
    xd = zeros(n_points)
    yd = zeros(n_points)

    for i in 1:n_points
        xs[i]=src[i][1]
        ys[i]=src[i][2]
        xd[i]=des[i][1]
        yd[i]=des[i][2]
    end

    #if homography can be exactly determined
    if n_points==4
        A = zeros(2*n_points,8)
        A[1:4,1]=xs
        A[1:4,2]=ys
        A[1:4,3]=1
        A[5:8,4]=xs
        A[5:8,5]=ys
        A[5:8,6]=1
        A[1:4,7]=-xs.*xd
        A[5:8,7]=-xs.*yd
        A[1:4,8]=-ys.*xd
        A[5:8,8]=-ys.*yd

        B = zeros(2*n_points)
        B[1:4]=xd
        B[5:8]=yd

        x = inv(A)*B

        push!(x,1)
        H = transpose(reshape(x, 3, 3))
        return H
    end

    A = zeros(2*n_points,9)
    A[1:n_points,1] = xs
    A[1:n_points,2] = ys
    A[1:n_points,3] = 1
    A[n_points+1:2*n_points,4]=xs
    A[n_points+1:2*n_points,5]=ys
    A[n_points+1:2*n_points,6]=1
    A[1:n_points,7] = -xs.*xd
    A[n_points+1:2*n_points,7]= -xs.*yd
    A[1:n_points,8] = -ys.*xd
    A[n_points+1:2*n_points,8]= -ys.*yd
    A[1:n_points,9] = -xd
    A[n_points+1:2*n_points,9]= -yd

    _, _, v = svd(A, thin=true)
    H=transpose(reshape(v[:,end]./v[end][end],3,3))

    return H
end