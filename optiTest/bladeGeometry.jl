### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 35e5a1ea-1aed-11ef-3356-8139a06354d9
begin
	import Pkg; Pkg.activate(".")
	using Plots, PlutoUI
	import FreeType, FileIO
end

# ╔═╡ a94613f8-147e-481d-91f3-4fe321e7cc99
airfoilCoords = [ 1.000000   .000000
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


# ╔═╡ c9bf9617-ff70-4be6-a791-368ba1ceb0b6
propgeom = [
    0.15   0.130   32.76
    0.20   0.149   37.19
    0.25   0.173   33.54
    0.30   0.189   29.25
    0.35   0.197   25.64
    0.40   0.201   22.54
    0.45   0.200   20.27
    0.50   0.194   18.46
    0.55   0.186   17.05
    0.60   0.174   15.97
    0.65   0.160   14.87
    0.70   0.145   14.09
    0.75   0.128   13.39
    0.80   0.112   12.84
    0.85   0.096   12.25
    0.90   0.081   11.37
    0.95   0.061   10.19
    1.00   0.041   8.99
    ]

# ╔═╡ bda62edf-313f-4d72-84e2-95bff2bc0a4b
begin
	Rtip = 0.127
	N = 2
	
	rs = propgeom[:, 1] * Rtip
	chords = propgeom[:, 2] * Rtip
	thetas = propgeom[:, 3] * -pi/180
end

# ╔═╡ 550b88a1-3746-424d-8a54-67fe84efc0f9
begin
	rotMat(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
	rot3Dyaxis(θ) = [cos(θ) 0 sin(θ);0 1 0; -sin(θ) 0 cos(θ)]
end

# ╔═╡ 497c3874-35bd-438a-b1c2-262b6358913a
begin
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
	    #AirfoilCoords *= [0 1; 1 0]
	    AirfoilCoords *= [1 0 0; 0 1 0]
	    AirfoilCoords[:,3] .= r
	    AirfoilCoords
	end
end

# ╔═╡ d932794c-d232-4958-8768-9b1cea786778
#rotate = @bind rotate Slider(1:.01:10, show_value=true)

# ╔═╡ a5921025-1a47-48a2-8abb-9089a55df95a
begin

	anim = @animate for rotate in 0:.01:10
		gr()
		#unicodeplots()
		theme(:solarized_light)
	    modAirfoilCoords = structAirfoilCoords(airfoilCoords, thetas, chords, rs)
	
	    modAirfoilCoords1 = similar(modAirfoilCoords)
	    modAirfoilCoords2 = similar(modAirfoilCoords); modAirfoilCoords3 = similar(modAirfoilCoords)
	
		ztheta = 2pi/N
	
		[modAirfoilCoords1[i,:,:] = modAirfoilCoords[i,:,:] *rot3Dyaxis(ztheta - rotate) for i in 1:length(rs)]
	    [modAirfoilCoords2[i,:,:] = modAirfoilCoords[i,:,:] *rot3Dyaxis(2 *ztheta - rotate) for i in 1:length(rs)]
	    #[modAirfoilCoords3[i,:,:] = modAirfoilCoords[i,:,:] *rot3Dyaxis(3 *ztheta - rotate) for i in 1:length(rs)]
	
	
	    p=plot(size=(1000,1000), xlims=(-0.15,0.15), ylims=(-0.15,0.15), zlims=(-0.15,0.15))
		
	    # [plot!(p, modAirfoilCoords1[i,:,1], modAirfoilCoords1[i,:,2], modAirfoilCoords1[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, grid=nothing, axis=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
	    # [plot!(p, modAirfoilCoords2[i,:,1], modAirfoilCoords2[i,:,2], modAirfoilCoords2[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, grid=nothing, axis=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
	    # [plot!(p, modAirfoilCoords3[i,:,1], modAirfoilCoords3[i,:,2], modAirfoilCoords3[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, grid=nothing, axis=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
	    [plot!(p, modAirfoilCoords1[i,:,1], modAirfoilCoords1[i,:,2], modAirfoilCoords1[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
	    [plot!(p, modAirfoilCoords2[i,:,1], modAirfoilCoords2[i,:,2], modAirfoilCoords2[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
	    #[plot!(p, modAirfoilCoords3[i,:,1], modAirfoilCoords3[i,:,2], modAirfoilCoords3[i,:,3], aspect_ratio=:equal, c=:green, label=nothing, xlabel="x", ylabel="y") for i in 1:length(rs)]
	end

	anim
end

# ╔═╡ 389c9af6-82b0-4634-bc64-65c9bae3acd1
gif(anim, "BladeSG6043.gif")

# ╔═╡ Cell order:
# ╠═35e5a1ea-1aed-11ef-3356-8139a06354d9
# ╟─a94613f8-147e-481d-91f3-4fe321e7cc99
# ╟─c9bf9617-ff70-4be6-a791-368ba1ceb0b6
# ╠═bda62edf-313f-4d72-84e2-95bff2bc0a4b
# ╠═550b88a1-3746-424d-8a54-67fe84efc0f9
# ╠═497c3874-35bd-438a-b1c2-262b6358913a
# ╠═d932794c-d232-4958-8768-9b1cea786778
# ╠═a5921025-1a47-48a2-8abb-9089a55df95a
# ╠═389c9af6-82b0-4634-bc64-65c9bae3acd1
