import Pkg; Pkg.activate(".")

using BenchmarkTools, Test
using Interpolations, NLsolve
using StructArrays

#################################################

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
	(alpha) -> linear_interpolation(airfoil.alphas, airfoil.CLs)( alphaShifter(alpha) ), (alpha) -> linear_interpolation(airfoil.alphas, airfoil.CDs)( alphaShifter(alpha) )
end

function alphaShifter(alpha)
	if -180 < alpha < 180
		return alpha
	elseif alpha > 180
		return alpha -180
	elseif alpha < -180
		return alpha +180
	end
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
	CL = section.CL(alpha *180/pi)
	CD = section.CD(alpha *180/pi)
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
calc_a(ϕ, F, σ′, Cn) = ( 1 - 4*F*sin(ϕ) /(σ′*Cn) )^-1
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

function solve_R(Rfunc::Function; x0=-0.01)
	sol = nlsolve((phi)->Rfunc(phi[1]), [x0], method=:broyden)
    #186.462 μs (1113 allocations: 245.34 KiB)
	sol.zero[1]
end

function integrate(y::AbstractVector)
    sum(y)
end

#################################################
#################rotor & op######################

airfoil = include("SG6043_360_Polar_NREL_Format.jl")
airfoildata = AirfoilData(airfoil[:,1], airfoil[:,2], airfoil[:,3])

interptest = interpCLCDstruct(airfoildata)
scatter(-180:180, d->interptest[1](d), size=(1000,500))
scatter(-180:180, d->interptest[2](d), size=(1000,500))

Rtip = 0.127
Rhub = 0.1*Rtip
N = 3
rotor = Rotor(Rtip,Rhub,N)

Vinf = 5.0
Omega = 5400*pi/30
rho = 1.225
op = OpCond(Vinf, Omega, rho)

#################################################
####################test0########################

begin
# 0.80   0.112   -12.84
	begin
	r = 0.8 *Rtip
	chord = 0.112 *Rtip
	theta = -12.84/180*pi
	sigmaPrime = N*chord/(2*pi*r)
	# R(ϕ) = R(ϕ, θ, F, σ′, Uinf, Ω, r, CL(), CD())
	R(phi)=R(phi, theta,1,sigmaPrime,Vinf,Omega,r,interpCLCDstruct(airfoildata)[1], interpCLCDstruct(airfoildata)[2])
	end

phiStar = solve_R(R)
alpha = theta - phiStar
CL = interpCLCDstruct(airfoildata)[1](alpha *180/pi)
CD = interpCLCDstruct(airfoildata)[2](alpha *180/pi)
Cn = calc_Cn(phiStar, CL, CD)
Ct = calc_Ct(phiStar, CL, CD)
a = calc_a(phiStar, 1, sigmaPrime, Cn)
aprime = calc_a′(phiStar, 1, sigmaPrime, Ct)

# 	calculate W
W = calc_W(Vinf, Omega, r, a, aprime)

# 	calculate T'
Tprime = calc_T(N, Cn, rho, W, chord, r)

# 	calculate Q'
Qprime = calc_Q(N, Ct, rho, W, chord, r)

alpha*180/pi, Tprime, Qprime, W
end

#################################################
####################test1########################

begin
r = propgeom[1,1] *Rtip
chord = propgeom[1,2] *Rtip
theta = propgeom[1,3] *pi/180

section = Section(r,chord,theta,airfoildata,rotor)

sol = Solve(rotor, section, op)
sol.alpha* 180/pi, sol.Tprime, sol.Qprime, sol.W
end

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
thetas = propgeom[:,3] *pi/180

sections = Section.(rs,chords,thetas,Ref(airfoildata),Ref(rotor))

sol = Solve.(Ref(rotor), sections, Ref(op))

sola = StructArray(sol)

sola.alpha[1] *180/pi
sola.Tprime[1] 
sola.Qprime[1]
sola.W[1]

sum(sola.Tprime) *(rs[2]-rs[1])
sum(sola.Qprime) *(rs[2]-rs[1])
sum(sola.Qprime) *(rs[2]-rs[1]) *op.Omega
sola.alpha .*180/pi

using Plots

scatter(propgeom[:,1], sola.Tprime)
scatter(rs, sola.Ct)
scatter(rs, sola.aPrime) 
