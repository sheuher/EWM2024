import Pkg; Pkg.activate(".")

using BenchmarkTools, Test
using Interpolations, NLsolve

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
	(alpha) -> linear_interpolation(airfoil.alphas, airfoil.CLs)(alpha), (alpha) -> linear_interpolation(airfoil.alphas, airfoil.CDs)(alpha)
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

	#   define residual function R(ϕ) = R(ϕ, θ, F, σ′, Uinf, Ω, r)
	RFunc(phi) = R(phi, section.theta, section.F(section.r), section.sigmaPrime, op.Vinf, op.Omega, section.r) 
    #   solve for R(ϕ*) = 0
	phi = solve_R(RFunc)

    # 	calculate Cn, Ct, a, a' as function of ϕ*
	alpha = section.theta - phi
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

	return Solve(RFunc, alpha, theta, phi, Cn, Ct, a, aPrime, W, Tprime, Qprime)
end


###############


interpCL(α) = linear_interpolation(airfoil[:,1], airfoil[:,2])(α)
interpCD(α) = linear_interpolation(airfoil[:,1], airfoil[:,3])(α)

calc_Cn(ϕ, CL, CD) = CL *cos(ϕ) - CD *sin(ϕ)
calc_Ct(ϕ, CL, CD) = - CL *sin(ϕ) + CD *cos(ϕ)
calc_a(ϕ, F, σ′, Cn) = ( 1 - 4*F*sin(ϕ) /(σ′*Cn) )^-1
calc_a′(ϕ, F, σ′, Ct) = ( 4 *F*sin(ϕ)*cos(ϕ) /(σ′*Ct) - 1 )^-1

calc_W(Uinf, Ω, r, a, a′) = sqrt(Uinf^2 *(1-a)^2 + (Ω*r*(1 +a′))^2)
calc_T(N, Cn, ρ, W, c, r) = N*Cn*0.5*ρ*W^2*c*r
calc_Q(N, Ct, ρ, W, c, r) = N*Ct*0.5*ρ*W^2*c*r^2


function R(ϕ, θ, F, σ′, Uinf, Ω, r)
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


function solve_R(Rfunc; x0=-0.01)
	sol = nlsolve((phi)->Rfunc(phi[1]), [x0], method=:broyden)
    #186.462 μs (1113 allocations: 245.34 KiB)
	print(sol)
	sol.zero[1]
end

####################


airfoil = include("naca4412.jl")
airfoil1 = AirfoilData(airfoil[:,1], airfoil[:,2], airfoil[:,3])


# geometry settings
begin
	Rtip = 0.127
	Rhub = 0.1 *Rtip
	N = 3
#
r = 0.95 *Rtip
chord = 0.061 *Rtip
theta = -10.19 *pi/180
sigmaPrime = N*chord/(2*pi*r)
end

# operating conditions
begin
	f = 1.0
	Vinf = 35.0
	Omega = 5400*pi/30
	rho = 1.225

	R((ϕ)) = R(ϕ, theta, f, sigmaPrime, Vinf, Omega, r) 
end

rotor1 = Rotor(Rtip, Rhub, N)
sectionInfo = Section(r,chord,theta,airfoil1,rotor1)
op = OpCond(Vinf, Omega, rho)
sol = Solve(rotor1, sectionInfo, op)
sol.alpha*180/pi, sol.Tprime, sol.Qprime, sol.W

begin
	# for each sections (r/R)
	# 	solve R(ϕ*) = 0
	phiStar = solve_R(R)
	
	# 	calculate Cn, Ct, a, a' as function of ϕ*
    alpha = theta - phiStar
	CL = interpCL(alpha*180/pi)
	CD = interpCD(alpha*180/pi)
	Cn = calc_Cn(phiStar, CL, CD)
	Ct = calc_Ct(phiStar, CL, CD)
	a = calc_a(phiStar, f, sigmaPrime, Cn)
	aprime = calc_a′(phiStar, f, sigmaPrime, Ct)
	
	# 	calculate W
	W = calc_W(Vinf, Omega, r, a, aprime)
	
	# 	calculate T'
	Tprime = calc_T(N, Cn, rho, W, chord, r)
	
	# 	calculate Q'
	Qprime = calc_Q(N, Ct, rho, W, chord, r)

	alpha*180/pi, Tprime, Qprime, W
end
	
# calculate T and Q by integrating through each section
#(Rtip -Rhub)*Tprime
#(Rtip -Rhub)*Qprime

# calculate power as Q Ω
#(Rtip -Rhub)*Qprime*Omega

let
	c = 1
	println("theta = $(theta *180/pi)")
	println("alpha = $(alpha *180/pi)")
	println("phiStar = $(phiStar *180/pi)")
	println("W = $(W)")
end