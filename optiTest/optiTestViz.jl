import Pkg; Pkg.activate(".")

using Plots
using DelimitedFiles

OptiRes = readdlm("BEMTOptiTestRUN.dat")
airfoilCoords = [
    1.000000   .000000
    .998105   .000656
    .992735   .002712
    .984387   .006072
    .973434   .010465
    .960071   .015523
    .944288   .020916
    .925966   .026547
    .905161   .032475
    .882072   .038683
    .856884   .045099
    .829794   .051648
    .801008   .058240
    .770741   .064777
    .739215   .071138
    .706663   .077173
    .673195   .082684
    .638892   .087606
    .603962   .091903
    .568537   .095505
    .532761   .098415
    .496847   .100592
    .460955   .102010
    .425276   .102690
    .389995   .102599
    .355272   .101752
    .321306   .100174
    .288269   .097880
    .256330   .094905
    .225671   .091274
    .196447   .087021
    .168817   .082199
    .142926   .076839
    .118895   .070995
    .096866   .064730
    .076936   .058094
    .059198   .051170
    .043757   .043981
    .030622   .036590
    .019828   .029159
    .011416   .021706
    .005277   .014355
    .001496   .007410
    .000024   .000942
    .000588  -.004677
    .003930  -.008595
    .010772  -.011195
    .020868  -.013164
    .034086  -.014302
    .050504  -.014598
    .070145  -.014215
    .092915  -.013293
    .118637  -.011970
    .147117  -.010343
    .178147  -.008477
    .211498  -.006446
    .246927  -.004303
    .284172  -.002104
    .322960   .000114
    .363008   .002319
    .404014   .004498
    .445705   .006689
    .487849   .008933
    .530290   .011167
    .572862   .013270
    .615252   .015037
    .657041   .016351
    .697849   .017212
    .737347   .017604
    .775201   .017514
    .811085   .016945
    .844678   .015916
    .875671   .014467
    .903773   .012654
    .928712   .010558
    .950241   .008276
    .968141   .005908
    .982146   .003613
    .992088   .001700
    .998026   .000453
    .999999   .000000]
begin    
    n = 20
    x = OptiRes
    rs = x[1:n]
    chords = x[n+1:2*n]
    thetas = x[2*n+1:3*n]

    Rhub = x[3*n+1]
    Rtip = x[3*n+2]

    Vinf = x[3*n+3]
    Omega = x[3*n+4]
    rho = x[3*n+5]
end

rotMat(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

plot(airfoilCoords[:,1], airfoilCoords[:,2], rs[1]*ones(35), size=(1000,500),  label=nothing, aspect_ratio=:equal)

function structAirfoilCoords(AirfoilCoords::AbstractMatrix, thetas::AbstractVector, chords::AbstractVector, rs::AbstractVector)
    @assert size(thetas) == size(chords) == size(rs)
    n = length(thetas)
    m = size(AirfoilCoords)[1] 
    res = Array{Float64}(undef, n, m, 3)
    [res[i, :,:] = structAirfoilCoords(AirfoilCoords, thetas[i], chords[i], rs[i]) for i in 1:n]
    res
end

function structAirfoilCoords(Coords::AbstractMatrix, theta::Float64, chord::Float64, r::Float64; cm::Float64=0.5)
    AirfoilCoords = copy(Coords)
    AirfoilCoords[:,1] .-= cm
    AirfoilCoords .*= chord
    AirfoilCoords *= rotMat(theta)
    AirfoilCoords *= [-1 0; 0 1]
    AirfoilCoords *= [1 0 0; 0 1 0]
    AirfoilCoords[:,3] .= r
    AirfoilCoords
end


let
coord1 = structAirfoilCoords(airfoilCoords, -0., 10., 2.)
coord2 = structAirfoilCoords(airfoilCoords, -15/180*pi, 10., 2.)

plot(coord1[:, 1], coord1[:,2], xlabel="x", ylabel="y", aspect_ratio=:equal, label="θ=0")
plot!(coord2[:, 1], coord2[:,2], xlabel="x", ylabel="y", aspect_ratio=:equal, label="θ=-15")
end

let
plotly()

modAirfoilCoords = structAirfoilCoords(airfoilCoords, thetas, chords, rs)
p=plot()
[plot!(p, modAirfoilCoords[i,:,1], modAirfoilCoords[i,:,2], modAirfoilCoords[i,:,3], aspect_ratio=:equal, label=nothing, grid=nothing, xlabel="x", ylabel="y") for i in 1:18]
p
end

let
    plotly()
    modAirfoilCoords = structAirfoilCoords(airfoilCoords, thetas, chords, rs)

    p=plot(size=(1000,1000), aspect_ratio=:equal)
    [plot!(p, modAirfoilCoords[i,:,1], modAirfoilCoords[i,:,2], modAirfoilCoords[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, xlabel="x", ylabel="y") for i in 1:18]
    #plot!(p, modAirfoilCoords[1,:,1], modAirfoilCoords[1,:,3])
    p
end

rot3Dyaxis(θ) = [cos(θ) 0 sin(θ);0 1 0; -sin(θ) 0 cos(θ)]
plotly()
let
    
    modAirfoilCoords = structAirfoilCoords(airfoilCoords, thetas, chords, rs)
	
    modAirfoilCoords1 = similar(modAirfoilCoords)
    modAirfoilCoords2 = similar(modAirfoilCoords); modAirfoilCoords3 = similar(modAirfoilCoords)

    ztheta = 2pi/3

    [modAirfoilCoords1[i,:,:] = modAirfoilCoords[i,:,:] *rot3Dyaxis(ztheta) for i in 1:length(rs)]
    [modAirfoilCoords2[i,:,:] = modAirfoilCoords[i,:,:] *rot3Dyaxis(2 *ztheta) for i in 1:length(rs)]
    [modAirfoilCoords3[i,:,:] = modAirfoilCoords[i,:,:] *rot3Dyaxis(3 *ztheta) for i in 1:length(rs)]


    p=plot(size=(1000,1000))
    
    # [plot!(p, modAirfoilCoords1[i,:,1], modAirfoilCoords1[i,:,2], modAirfoilCoords1[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, grid=nothing, axis=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
    # [plot!(p, modAirfoilCoords2[i,:,1], modAirfoilCoords2[i,:,2], modAirfoilCoords2[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, grid=nothing, axis=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
    # [plot!(p, modAirfoilCoords3[i,:,1], modAirfoilCoords3[i,:,2], modAirfoilCoords3[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, grid=nothing, axis=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
    [plot!(p, modAirfoilCoords1[i,:,1], modAirfoilCoords1[i,:,2], modAirfoilCoords1[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
    [plot!(p, modAirfoilCoords2[i,:,1], modAirfoilCoords2[i,:,2], modAirfoilCoords2[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
    [plot!(p, modAirfoilCoords3[i,:,1], modAirfoilCoords3[i,:,2], modAirfoilCoords3[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
    #savefig(p, "propeller0.html")
    p
end

plot(thetas*180/pi)
plot(chords)
Omega
