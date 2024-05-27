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
		(alpha) -> begin
			linear_interpolation(airfoil.alphas, airfoil.CLs)(alphaShifter(alpha))
		end, (alpha) -> begin
			linear_interpolation(airfoil.alphas, airfoil.CDs)(alphaShifter(alpha))
		end
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
	Section(r,chord,theta,airfoil::AirfoilData,rotor::Rotor) = Section( r,chord,theta,
	                                                                    interpCLCDstruct(airfoil)..., 
	                                                                    (r)->1, 
	                                                                    rotor.N*chord/(2*pi*r) )
	
	
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
	    #   solve for R(ϕ*) = 0
		phi = solve_R(RFunc)
	
	    # 	calculate Cn, Ct, a, a' as function of ϕ*
		alpha = section.theta - phi
		alpha = alphaShifter(alpha)
		CL = section.CL(alpha)
		CD = section.CD(alpha)
		Cn = calc_Cn(phi, CL, CD)
		Ct = calc_Ct(phi, CL, CD)
		a = calc_a(phi, section.F(section.r), section.sigmaPrime, Cn)
		aPrime = calc_a′(phi, section.F(section.r), section.sigmaPrime, Ct)
		
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
	calc_Ct(ϕ, CL, CD) = - CL *sin(ϕ) + CD *cos(ϕ)
	calc_a(ϕ, F, σ′, Cn) = ( 1 - 4*F*(sin(ϕ))^2 /(σ′*Cn) )^-1
	calc_a′(ϕ, F, σ′, Ct) = ( 4 *F*sin(ϕ)*cos(ϕ) /(σ′*Ct) - 1 )^-1
	
	calc_W(Uinf, Ω, r, a, a′) = sqrt(Uinf^2 *(1 -a)^2 + (Ω*r*(1 +a′))^2)
	calc_T(N, Cn, ρ, W, c, r) = N*Cn*0.5*ρ*W^2*c*r
	calc_Q(N, Ct, ρ, W, c, r) = N*Ct*0.5*ρ*W^2*c*r^2
	
	
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
		a = calc_a(ϕ, F, σ′, Cn)
	
		# calculate a' (according to Nings2020 derivation)
		a′ = calc_a′(ϕ, F, σ′, Ct)
		
		# return Residual
		return sin(ϕ)/(a -1) - Uinf*cos(ϕ) /(Ω*r*(1-a′)) 
	end
	
	function solve_R(Rfunc::Function)
		# using Roots's Brent-Dekker method
		if (Rfunc(-pi/2) > 0) & (Rfunc(-pi/2)*Rfunc(-1e-6) < 0)
			fzero(Rfunc, (-pi/2, -1e-6), Roots.Brent())
		else
		#	fzero(Rfunc, (-pi/2, -pi), Roots.Brent())
			999
		end
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
	p1 = scatter(-10pi:0.01:10pi, d->interptest[1](d), size=(1000,500))
	p2 = scatter(-2pi:0.01:2pi, d->interptest[2](d), size=(1000,500))

	plot(p1,p2,layout=(2,1))
end

# ╔═╡ b38ebafd-c6e2-4e18-8d63-8b943d7d439b
Vinf = @bind Vinf Slider(1:35, show_value=true)

# ╔═╡ bef3ad20-9a9f-4cf2-bf3a-5ffa1b8dc785
begin
		
	#################################################
	#################rotor & op######################
		
		Rtip = 0.127
		Rhub = 0.1*Rtip
		N = 3
		rotor = Rotor(Rtip,Rhub,N)
		
		#Vinf = 5.0
		Omega = 5400*pi/30
		rho = 1.225
		op = OpCond(Vinf, Omega, rho)
		
end

# ╔═╡ bb229950-7eaf-4c4c-aa61-9aee13fc1c09
begin
	#################################################
	#####################test2#######################
	
	
	#   r/R, c/R, theta[degree]
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
	
	rs = propgeom[:,1] *Rtip
	chords = propgeom[:,2] *Rtip
	thetas = - propgeom[:,3] *pi/180
	
	sections = Section.(rs,chords,thetas,Ref(airfoildata),Ref(rotor))
	
	sol = Solve.(Ref(rotor), sections, Ref(op))
	
	sola = StructArray(sol)
	
	#sola.alpha[1] *180/pi
	#sola.Tprime[1] 
	#sola.Qprime[1]
	#sola.W[1]
	
	println("Thrust = $(sum(sola.Tprime) *(rs[2]-rs[1]))")
	println("Qprime = $(sum(sola.Qprime) *(rs[2]-rs[1]))")
	println("Power = $(sum(sola.Qprime) *(rs[2]-rs[1]) *op.Omega)")
	#sola.alpha .*180/pi
	
	
	p5 = scatter(propgeom[:,1], sola.Qprime, ylims=(-0.8,2), ms=2, label="Q'")
	p6 = scatter(rs, sola.alpha *180/pi, ms=1, label="alpha")
	p7 = scatter(rs, sola.theta *180/pi, ms=1, label="theta")
	p4 = scatter(rs, sola.phi *180/pi, ms=1, label="phi") 
	plot(p7,p4,p6, p5,layout=(4,1))

end

# ╔═╡ e446c57f-e1e4-4616-a706-85ac14bbcd86
sola.RFunction[1](-pi/2)

# ╔═╡ 87346d7f-8a7d-42b0-b95b-fa6d6bf52b5a
1-0.3

# ╔═╡ Cell order:
# ╠═d9b8d3ee-1b2b-11ef-3574-1991aad1e852
# ╠═ecdc3a49-fc13-4f23-b145-edb406629a18
# ╠═67cb8ae8-4d67-497f-ac9a-ad316336ca8e
# ╠═ef399a25-87b2-457c-83c3-7cb757ed74f8
# ╠═bef3ad20-9a9f-4cf2-bf3a-5ffa1b8dc785
# ╠═b38ebafd-c6e2-4e18-8d63-8b943d7d439b
# ╠═bb229950-7eaf-4c4c-aa61-9aee13fc1c09
# ╠═e446c57f-e1e4-4616-a706-85ac14bbcd86
# ╠═87346d7f-8a7d-42b0-b95b-fa6d6bf52b5a
