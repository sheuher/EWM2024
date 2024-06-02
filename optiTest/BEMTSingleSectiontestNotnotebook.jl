begin
	import Pkg; Pkg.activate(".")
	
	using Interpolations
	using StructArrays
	using Plots
	using DelimitedFiles
	using NonlinearSolve, Roots, NLsolve
end

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
end

begin
	airfoil = readdlm("SG6043_360_Polar_NREL_Format.txt")
	airfoildata = AirfoilData(airfoil[:,1]*pi/180, airfoil[:,2], airfoil[:,3])
end


begin
	interptest = interpCLCDstruct(airfoildata)
	p1 = scatter(-2pi:0.01:2pi, d->interptest[1](d), size=(1000,500), ms=1)
	p2 = scatter(-2pi:0.01:2pi, d->interptest[2](d), size=(1000,500), ms=1)

	plot(p1,p2,layout=(2,1))
end




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

function turbine(Ω, Vinf; Ncut=102)
	Rtip = 0.55
	Rhub = 0.135
	N = 3
	rotor = Rotor(Rtip,Rhub,N)
	
	Vinf = Vinf
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

turbineZero(omega, v) = turbine(omega, v; Ncut=34).M - 3.7 

# solve for M(omega) - 3.7 = 0
nlsolve((omega) -> turbineZero(omega[1], 5), [250.]).zero

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


begin
	omegas = 12:500
	Vinfty = 10
	resTot = turbine.( omegas, Vinfty; Ncut=104 )
end

let sola=resTot[1].sol
	p5 = scatter(rs, sola.Qprime,  ms=2, label="Q'")
	p6 = scatter(rs, sola.alpha *180/pi, ms=1, label="alpha")
	p7 = scatter(rs, sola.theta *180/pi, ms=1, label="theta")
	p4 = scatter(rs, sola.phi *180/pi, ms=1, label="phi") 
	plot(p7,p4,p6, p5,layout=(4,1), size=(750,750))
end



resTotS = StructArray(resTot)

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

# Thrust = ~ 10.9
# Torque = ~ 2.8
# Power = ~ 3.4 
# with 5 elements, for v=10, Omega=12
