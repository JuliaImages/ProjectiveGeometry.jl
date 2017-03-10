module CameraGeometryTest

# write your own tests here
using FactCheck, CameraGeometry, Base.Test, Images

include("homography.jl")

isinteractive() || FactCheck.exitstatus()

end
