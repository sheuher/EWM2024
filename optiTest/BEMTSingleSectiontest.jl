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

# ╔═╡ d9b8d3ee-1b2b-11ef-3574-1991aad1e852
begin
	import Pkg; Pkg.activate(".")
	
	using PlutoUI
	using Interpolations, NLsolve
	using StructArrays
	using Plots
	using DelimitedFiles
	using NonlinearSolve, Roots
end

# ╔═╡ ecdc3a49-fc13-4f23-b145-edb406629a18
begin
	
	struct AirfoilData
	    alphas :: AbstractVector
	    CLs :: AbstractVector
	    CDs :: AbstractVector
	end
	
	struct Rotor
		Rtip :: Float64
		Rhub :: Float64
		N :: Int64
	end
	
	function interpCLCDstruct(airfoil::AirfoilData)
		interpStruct(airfoil.alphas, airfoil.CLs), interpStruct(airfoil.alphas, airfoil.CDs)
	end

	function interpStruct(alphas::AbstractVector, y::AbstractVector)
		# linear interpolation with periodic phase shift correction 
		(x) -> linear_interpolation(alphas,y)(alphaShifter(x))
	end
	
	function alphaShifter(alpha)
		asin(sin(alpha))
	end
		
	struct Section
		r :: Float64
		chord :: Float64
		theta :: Float64
		CL :: Function
		CD :: Function
		F :: Function
		sigmaPrime :: Float64
	end
	Section(r,chord,theta,airfoil::AirfoilData,rotor::Rotor) = Section( r,chord,theta,interpCLCDstruct(airfoil)..., (r)->1, rotor.N*chord/(2*pi*r) )
	
	
	struct OpCond
		Vinf :: Float64
		Omega :: Float64
		rho :: Float64
	end
	
	mutable struct Solve
		RFunction :: Function
		alpha :: Float64
		theta :: Float64
		phi :: Float64
		Cn :: Float64
		Ct :: Float64
		a :: Float64
		aPrime :: Float64
		W :: Float64
		Tprime :: Float64
		Qprime :: Float64
	end
	
	function Solve(rotor::Rotor, section::Section, op::OpCond)
	
		#   define residual function R(ϕ) = R(ϕ, θ, F, σ′, Uinf, Ω, r, CL(), CD())
		RFunc(phi) = R(phi, section.theta, section.F(section.r), section.sigmaPrime, op.Vinf, op.Omega, section.r, section.CL, section.CD) 

		#   define residual function RPB(ϕ) = RPB(ϕ, θ, F, σ′, Uinf, Ω, r, CL(), CD()) for propeller brake zone
		RPBFunc(phi) = RPB(phi, section.theta, section.F(section.r), section.sigmaPrime, op.Vinf, op.Omega, section.r, section.CL, section.CD)
		
	    #   solve for R(ϕ*) = 0
		phi = solve_BEM(RFunc, RPBFunc)
	
	    # 	calculate Cn, Ct, a, a' as function of ϕ*
		alpha = section.theta - phi
		alpha = alphaShifter(alpha)
		CL = section.CL(alpha)
		CD = section.CD(alpha)
		Cn = calc_Cn(phi, CL, CD)
		Ct = calc_Ct(phi, CL, CD)
		
		kappa = calc_kappa(phi, section.F(section.r), section.sigmaPrime, Cn)
		a = calc_a(kappa, section.F(section.r))
		kappaPrime = calc_kappaPrime(phi, section.F(section.r), section.sigmaPrime, Ct)
		aPrime = calc_aPrime(kappaPrime)
		
		# 	calculate W
		W = calc_W(op.Vinf, op.Omega, section.r, a, aPrime)
		
		# 	calculate T'
		Tprime = calc_T(rotor.N, Cn, op.rho, W, section.chord, section.r)
		
		# 	calculate Q'
		Qprime = calc_Q(rotor.N, Ct, op.rho, W, section.chord, section.r)
	
		Solve(RFunc, alpha, section.theta, phi, Cn, Ct, a, aPrime, W, Tprime, Qprime)
	end
	
	#################################################
	
	calc_Cn(ϕ, CL, CD) = CL *cos(ϕ) - CD *sin(ϕ)
	calc_Ct(ϕ, CL, CD) = -CL *sin(ϕ) - CD *cos(ϕ)

	calc_kappa(ϕ, F, σ′, Cn) = (σ′*Cn) / (4*F*(sin(ϕ))^2)
	calc_kappaPrime(ϕ, F, σ′, Ct) = (-σ′*Ct) / (4 *F*sin(ϕ)*cos(ϕ))
	calc_a(κ) = κ / (1 +κ)
	calc_aPrime(κ′) = κ′ / (1 -κ′)

	function calc_a(κ, F)
		if κ < 2/3
			# momentum for windmill
			return κ / (1 +κ) #calc_a(κ)
		else
			# empirical with Buhl derivation 
			γ₁ = 2*F*κ - (10/9 -F)
			γ₂ = 2*F*κ - F*(4/3 -F)
			γ₃ = 2*F*κ - (25/9 -2*F)
			println(2)
			return (γ₁ -√(γ₂)) / γ₃
		end
	end
	
	calc_W(Uinf, Ω, r, a, a′) = sqrt(Uinf^2 *(1 -a)^2 + (Ω*r*(1 +a′))^2)
	calc_T(N, Cn, ρ, W, c, r) = N*Cn*0.5*ρ*W^2*c # N direction
	calc_Q(N, Ct, ρ, W, c, r) = N*Ct*0.5*ρ*W^2*c*r # T direction

	function calc_F()
	end
	
	
	function R(ϕ, θ, F, σ′, Uinf, Ω, r, interpCL::Function, interpCD::Function)
		# calculate angle of attack
		α = θ - ϕ
	
		# calculate CL = f(α)
		CL = interpCL(α)
	
		# calculate CD = f(α)
		CD = interpCD(α)
	
		# calculate Cn as ϕ
		Cn = calc_Cn(ϕ, CL, CD)
	
		# calculate Ct as ϕ
		Ct = calc_Ct(ϕ, CL, CD) 
	
		# calculate a (according to Nings2020 derivation)
		κ = calc_kappa(ϕ, F, σ′, Cn)
		a = calc_a(κ, F)
	
		# calculate a' (according to Nings2020 derivation)
		κ′ = calc_kappaPrime(ϕ, F, σ′, Ct)
		a′ = calc_aPrime(κ′)
		
		# return Residual
		return sin(ϕ)/(a -1) - Uinf*cos(ϕ) /(Ω*r*(1 +a′)) 
	end

	function RPB(ϕ, θ, F, σ′, Uinf, Ω, r, interpCL::Function, interpCD::Function)
		# calculate angle of attack
		α = θ - ϕ
	
		# calculate CL = f(α)
		CL = interpCL(α)
	
		# calculate CD = f(α)
		CD = interpCD(α)
	
		# calculate Cn as ϕ
		Cn = calc_Cn(ϕ, CL, CD)
	
		# calculate Ct as ϕ
		Ct = calc_Ct(ϕ, CL, CD) 
	
		# calculate a (according to Nings2020 derivation)
		κ = calc_kappa(ϕ, F, σ′, Cn)
		a = calc_a(κ, F)
	
		# calculate a' (according to Nings2020 derivation)
		κ′ = calc_kappaPrime(ϕ, F, σ′, Ct)
		a′ = calc_aPrime(κ′)
		
		# return Residual
		return sin(ϕ)/(κ -1) - Uinf*cos(ϕ)*(1 -κ′) /(Ω*r) 
	end
	
	function solve_R(Rfunc::Function)
		# using Roots's Brent-Dekker method
		fzero(Rfunc, (-pi/2, -1e-6), Roots.Brent())
	end

	function solve_R(Rfunc::Function, lb, ub)
		# using Roots's Brent-Dekker method
		fzero(Rfunc, (lb, ub), Roots.Brent())
	end

	function solve_BEM(Rfunc::Function, RPBfunc::Function; ϵ=-1e-6)
		if Rfunc(-π/2) > 0
			return solve_R(Rfunc,-π/2,ϵ)
		elseif (RPBfunc(-π/4) < 0) && (RPBfunc(ϵ) > 0)
			println(4)
			return solve_R(RPBfunc,-π/4,-ϵ)
		else
			println(5)
			return solve_R(Rfunc,π/2,π)
		end
	end

	function trapz(x, y)

	    integral = 0.0
	    for i = 1:length(x)-1
	        integral += (x[i+1]-x[i])*0.5*(y[i] + y[i+1])
	    end
	    return integral
	end
end

# ╔═╡ 67cb8ae8-4d67-497f-ac9a-ad316336ca8e
begin
	airfoil = readdlm("SG6043_360_Polar_NREL_Format.txt")
	airfoildata = AirfoilData(airfoil[:,1]*pi/180, airfoil[:,2], airfoil[:,3])
end

# ╔═╡ ef399a25-87b2-457c-83c3-7cb757ed74f8
begin
	interptest = interpCLCDstruct(airfoildata)
	p1 = scatter(-2pi:0.01:2pi, d->interptest[1](d), size=(1000,500), ms=1)
	p2 = scatter(-2pi:0.01:2pi, d->interptest[2](d), size=(1000,500), ms=1)

	plot(p1,p2,layout=(2,1))
end

# ╔═╡ bef3ad20-9a9f-4cf2-bf3a-5ffa1b8dc785
begin
		
	#################################################
	#################rotor & op######################
		
		Rtip = 0.55
		Rhub = 0.135
		N = 3
		rotor = Rotor(Rtip,Rhub,N)
		
		Vinf = 10
		Omega = 120/60*2*pi
		rho = 1.225
		op = OpCond(Vinf, Omega, rho)
	op.Omega
		
end

# ╔═╡ b38ebafd-c6e2-4e18-8d63-8b943d7d439b
#Vinf = @bind Vinf Slider(1:35, show_value=true)

# ╔═╡ bb229950-7eaf-4c4c-aa61-9aee13fc1c09
begin
	#################################################
	#####################test2#######################
	
	
	#   r/R, c/R, theta[degree]
	propgeom = [
	    0.135/0.55   0.20   -30
	    1.00   0.15   -9
	    ]
	
	rs = LinRange( propgeom[1,1], propgeom[2,1], 102 )[2:end-1] *Rtip
	chords = LinRange( propgeom[1,2], propgeom[2,2], 102 )[2:end-1] *Rtip
	thetas = LinRange( propgeom[1,3], propgeom[2,3], 102 )[2:end-1] *pi/180
end

# ╔═╡ 06546b4a-b706-425e-918a-de0f29a0bc29
begin
	sections = Section.(rs,chords,thetas,Ref(airfoildata),Ref(rotor))
	
	sol = Solve.(Ref(rotor), sections, Ref(op))
	
	sola = StructArray(sol)
	
	#sola.alpha[1] *180/pi
	#sola.Tprime[1] 
	#sola.Qprime[1]
	#sola.W[1]

	KP = -sections[1].sigmaPrime * sola[1].Ct / (4*cos(sola[1].phi)*sin(sola[1].phi))
	calc_aPrime(KP)
end

# ╔═╡ 19e07920-91a6-4150-9104-cd91b9fbd24b
begin
	function ThrustTotal(sols, sections, rotor, op)
		# numerical integration for thrust for linearly spaced radial section
		rs = [section.r for section in sections]
		Tfull = [0.; sols.Tprime; 0]
		rfull = [rotor.Rhub; rs; rotor.Rtip]
		sum(Tfull) * (rfull[3] - rfull[2])
	end
	function TorqueTotal(sols, sections, rotor, op)
		# numerical integration for thrust for linearly spaced radial section
		rs = [section.r for section in sections]
		Qfull = [0.; sols.Qprime; 0]
		rfull = [rotor.Rhub; rs; rotor.Rtip]
		sum(Qfull) * (rfull[3] - rfull[2])
	end
	function PowerTotal(Torque, Omega)
		Torque *Omega
	end
end

# ╔═╡ ed5d4dd3-af9d-4722-8bc2-54730dd9df48
begin
	Tfull = [0.; sola.Tprime; 0]
	Qfull = [0.; sola.Qprime; 0]
	
	rfull = [rotor.Rhub; rs; rotor.Rtip]
	Thrust = ThrustTotal(sola, sections, rotor, op)#sum(Tfull) * (rfull[3] - rfull[2])
	Torque = TorqueTotal(sola, sections, rotor, op)#sum(Qfull) * (rfull[3] - rfull[2])

	
	Power = PowerTotal(Torque, op.Omega)#Torque * op.Omega
	
	println("Thrust = $(Thrust)")
	println("Torque = $(Torque)")
	println("Power = $(Power)")
	#sola.alpha .*180/pi
	
	
	p5 = scatter(rs, sola.Qprime,  ms=2, label="Q'")
	p6 = scatter(rs, sola.alpha *180/pi, ms=1, label="alpha")
	p7 = scatter(rs, sola.theta *180/pi, ms=1, label="theta")
	p4 = scatter(rs, sola.phi *180/pi, ms=1, label="phi") 
	plot(p7,p4,p6, p5,layout=(4,1), size=(750,750))

end

# ╔═╡ cb0c4980-1203-4328-a4a6-c657816ba785
function turbine(Ω, V; Ncut=102)
	Rtip = 0.55
	Rhub = 0.135
	N = 3
	rotor = Rotor(Rtip,Rhub,N)
	
	Vinf = V
	Omega = Ω/60*2*pi
	rho = 1.225
	OP = OpCond(Vinf, Omega, rho)
	propgeom = [
	    0.135/0.55   0.20   -30
	    1.00   0.15   -9
	    ]
	
	rs = LinRange( propgeom[1,1], propgeom[2,1], Ncut )[2:end-1] *Rtip
	chords = LinRange( propgeom[1,2], propgeom[2,2], Ncut )[2:end-1] *Rtip
	thetas = LinRange( propgeom[1,3], propgeom[2,3], Ncut )[2:end-1] *pi/180

	sections = Section.(rs,chords,thetas,Ref(airfoildata),Ref(rotor))
	
	sol = Solve.(Ref(rotor), sections, Ref(OP))
	
	sola = StructArray(sol)

	T = ThrustTotal(sola, sections, rotor, OP)
	M = TorqueTotal(sola, sections, rotor, OP)
	P = PowerTotal(M, OP.Omega)

	(sol = sola, T=T, M=M, P=P)
end

# ╔═╡ 2ea86e2b-63fd-4c9c-926c-3fd30a5674ba
turbineZero(omega, v) = turbine(omega, v; Ncut=34).M - 3.7 

# ╔═╡ 15a51040-b278-4f18-aa76-f01bbd11eba4
vvv = @bind vvv Slider(1.:0.2:18., show_value=true)

# ╔═╡ 9eeeaab2-33c7-4339-b179-f079ab144574
plot(0:20:1500, turbineZero.(0:20:1500, vvv))

# ╔═╡ a004b484-a861-45ef-9d93-d60db9150478
# solve for M(omega) - 3.7 = 0
nlsolve((omega) -> turbineZero(omega[1], 11), [250.]).zero

# ╔═╡ 67edacb0-528c-4066-a2e0-6a88a17cee7d
let
	OME=let
	Vinfs = 8.6:0.5:18
	map(Vinfs) do Vinf
		nlsolve((omega) -> turbineZero(omega[1], Vinf), [1000.]).zero[1]
	end
end

#surface(12:500, 3:15, turbineZero)

po=[turbine(OME[i], (8.6:0.5:18)[i]; Ncut=34).T for i in 1:length(OME)]

let
p1=plot(8.6:0.5:18, po, xlabel="Vinf", ylabel="Power")
p2=plot(8.6:0.5:18, OME, xlabel="Vinf", ylabel="Omega")
plot(p1,p2,layout=(2,1), size=(1000,1000), legend=nothing)
end

end

# ╔═╡ e131f2aa-510c-4d2e-b709-cca96fcab395
sqrt((258/60*2pi * 0.55)^2 + (10)^2) # W ≈ 17.91

# ╔═╡ 39d504bc-7ccd-412a-92da-d27dd4f48198
begin
	omegas = 12:500
	Vinfty = 10
	resTot = turbine.( omegas, Vinfty; Ncut=104 )
end

# ╔═╡ a374be1c-4630-4530-a055-d13e6228391e
let sola=resTot[1].sol
	p5 = scatter(rs, sola.Qprime,  ms=2, label="Q'")
	p6 = scatter(rs, sola.alpha *180/pi, ms=1, label="alpha")
	p7 = scatter(rs, sola.theta *180/pi, ms=1, label="theta")
	p4 = scatter(rs, sola.phi *180/pi, ms=1, label="phi") 
	plot(p7,p4,p6, p5,layout=(4,1), size=(750,750))
end

# ╔═╡ bf8e65aa-f629-4be7-9948-4ccf56050a88


# ╔═╡ 0670a6e0-ae8f-4b82-a9b6-c655e8b5d7f8


# ╔═╡ 61452a85-cb0d-4618-bf0f-0099209eb2b4
resTotS = StructArray(resTot)

# ╔═╡ 6dbf4c10-cf7f-4a8f-bd3c-284b2cc478e4
let
	p1 = plot(omegas, resTotS.M, xticks=0:20:500, yticks=minimum(resTotS.M):0.2:maximum(resTotS.M)+5, 
	xlabel="Ω [Hz]", ylabel="Torque [Nm]", label=nothing)
	p2 = plot(omegas, resTotS.T, xticks=0:20:500, yticks=minimum(resTotS.T):2:maximum(resTotS.T)+5,
		xlabel="Ω [Hz]", ylabel="Thrust [N]", label=nothing)
	p3 = plot(omegas, resTotS.P, xticks=0:20:500, yticks=minimum(resTotS.P):20:maximum(resTotS.P)+20,
		xlabel="Ω [Hz]", ylabel="Power [W]", label=nothing)
	plot(p1,p2,p3,layout=(3,1),size=(1000,1000))
	title!("V=$(Vinfty) m/s")
end

# ╔═╡ a31de0a2-62d9-47a1-ab1a-3087c6d0d166
# Thrust = ~ 10.9
# Torque = ~ 2.8
# Power = ~ 3.4 
# with 5 elements, for v=10

# ╔═╡ 5a9031f8-a279-4cb1-83c1-0a431c0c603a
#sola.aPrime |> plot # this could be the main wrong point, a is totally correct

# ╔═╡ b5f68bab-3be5-4884-9d6c-dca0c6277cc0
#sola.phi*180/pi |> plot# this could be the main wrong point

# ╔═╡ 5a6ceeb4-19d5-4d56-ba9f-53117cd00425
i = @bind i Slider(1:18, show_value=true)

# ╔═╡ e446c57f-e1e4-4616-a706-85ac14bbcd86
#plot(-pi/2:0.001:pi/2, sola.RFunction[i], ylims=(-5,5))

# ╔═╡ Cell order:
# ╠═d9b8d3ee-1b2b-11ef-3574-1991aad1e852
# ╠═ecdc3a49-fc13-4f23-b145-edb406629a18
# ╠═67cb8ae8-4d67-497f-ac9a-ad316336ca8e
# ╠═ef399a25-87b2-457c-83c3-7cb757ed74f8
# ╠═bef3ad20-9a9f-4cf2-bf3a-5ffa1b8dc785
# ╠═b38ebafd-c6e2-4e18-8d63-8b943d7d439b
# ╠═bb229950-7eaf-4c4c-aa61-9aee13fc1c09
# ╠═06546b4a-b706-425e-918a-de0f29a0bc29
# ╠═ed5d4dd3-af9d-4722-8bc2-54730dd9df48
# ╠═19e07920-91a6-4150-9104-cd91b9fbd24b
# ╠═cb0c4980-1203-4328-a4a6-c657816ba785
# ╠═2ea86e2b-63fd-4c9c-926c-3fd30a5674ba
# ╠═15a51040-b278-4f18-aa76-f01bbd11eba4
# ╠═9eeeaab2-33c7-4339-b179-f079ab144574
# ╠═a004b484-a861-45ef-9d93-d60db9150478
# ╠═67edacb0-528c-4066-a2e0-6a88a17cee7d
# ╟─e131f2aa-510c-4d2e-b709-cca96fcab395
# ╠═39d504bc-7ccd-412a-92da-d27dd4f48198
# ╠═a374be1c-4630-4530-a055-d13e6228391e
# ╠═bf8e65aa-f629-4be7-9948-4ccf56050a88
# ╠═0670a6e0-ae8f-4b82-a9b6-c655e8b5d7f8
# ╠═61452a85-cb0d-4618-bf0f-0099209eb2b4
# ╠═6dbf4c10-cf7f-4a8f-bd3c-284b2cc478e4
# ╠═a31de0a2-62d9-47a1-ab1a-3087c6d0d166
# ╠═5a9031f8-a279-4cb1-83c1-0a431c0c603a
# ╠═b5f68bab-3be5-4884-9d6c-dca0c6277cc0
# ╠═5a6ceeb4-19d5-4d56-ba9f-53117cd00425
# ╠═e446c57f-e1e4-4616-a706-85ac14bbcd86
