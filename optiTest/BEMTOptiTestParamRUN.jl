begin
	import Pkg; Pkg.activate(".")
	
	using Interpolations
	using StructArrays
	using Plots
	using DelimitedFiles
	#using NonlinearSolve, Roots, NLsolve
    using Ipopt
    using ForwardDiff
	using FLOWMath#:brent
end

begin
	struct AirfoilData
	    alphas :: AbstractVector
	    CLs :: AbstractVector
	    CDs :: AbstractVector
	end
	
	struct Rotor
		Rtip 
		Rhub 
		N 
	end
	
	function interpCLCDstruct(airfoil::AirfoilData)
		interpStruct(airfoil.alphas, airfoil.CLs), interpStruct(airfoil.alphas, airfoil.CDs)
	end

	function interpStruct(alphas, y)
		# linear interpolation with periodic phase shift correction 
		(x) -> linear_interpolation(alphas,y)(alphaShifter(x))
	end
	
	function alphaShifter(alpha)
		asin(sin(alpha))
	end
		
	struct Section
		r 
		chord 
		theta 
		CL :: Function
		CD :: Function
		F :: Function
		sigmaPrime 
	end
	Section(r,chord,theta,airfoil::AirfoilData,rotor::Rotor) = Section( r,chord,theta,interpCLCDstruct(airfoil)..., (r)->1, rotor.N*chord/(2*pi*r) )
	
	
	struct OpCond
		Vinf 
		Omega 
		rho 
	end
	
	mutable struct Solve
		RFunction :: Function
		alpha 
		theta 
		phi 
		Cn 
		Ct 
		a 
		aPrime 
		W 
		Tprime 
		Qprime 
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
		# tip loss
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
		#fzero(Rfunc, (lb, ub), Roots.Brent())

		# a better Brent has to be used here
		sol, _ = brent(Rfunc, lb, ub)
		sol
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
    nothing
end

#==========================================================================#

begin
	airfoil = readdlm("SG6043_360_Polar_NREL_Format.dat", skipstart=3)
	airfoildata = AirfoilData(airfoil[:,1], airfoil[:,2], airfoil[:,3])

	Rtip = 0.55
	Rhub = 0.135
	N = 3
	Vinf = 10.
	Nparam = 5
	Ncut = 50
	Omega = 500.0/60*2*pi
	rho = 1.225
	propgeom = [0.135/0.55   0.20   -30
				1.00         0.15   -9]

	#r = LinRange( Rhub/Rtip, 1., Ncut+2 )[2:end-1] *Rtip
	#chord = LinRange( propgeom[1,2], propgeom[2,2], Ncut+2 )[2:end-1] *Rtip
	#theta = LinRange( propgeom[1,3], propgeom[2,3], Ncut+2 )[2:end-1] *pi/180

	rP = LinRange( (1 -Rhub/Rtip)/(Ncut+1) +Rhub/Rtip, 1 -(1 -Rhub/Rtip)/(Ncut+1), Nparam ) *Rtip
	chordP = LinRange( propgeom[1,2], propgeom[2,2], Nparam ) *Rtip
	thetaP = LinRange( propgeom[1,3], propgeom[2,3], Nparam ) *pi/180

	r = LinRange( Rhub/Rtip, 1., Ncut+2 )[2:end-1] *Rtip
	chord = akima(rP, chordP, r)
	theta = akima(rP, thetaP, r)

	x = [rP; chordP; thetaP; Rhub; Rtip; Vinf; Omega; rho]
	y = zeros(2)

	function BEMTwrapper(x)
		n = Nparam
		rP = x[1:n]
		chordP = x[n+1:2*n]
		thetaP = x[2*n+1:3*n]

		Rhub = x[3*n+1]
		Rtip = x[3*n+2]

		Vinf = x[3*n+3]
		Omega = x[3*n+4]
		rho = x[3*n+5]

		r = LinRange( Rhub/Rtip, 1., Ncut+2 )[2:end-1] *Rtip
		chord = akima(rP, chordP, r)
		theta = akima(rP, thetaP, r)

		rotor = Rotor(Rtip, Rhub, N)
		sections = Section.(r, chord, theta, interpCLCDstruct(airfoildata)..., (r)->1, N*chord ./(2*pi*r))
		op = OpCond(Vinf, Omega, rho)

		sols = Solve.(Ref(rotor), sections, Ref(op))
		outputs = StructArray(sols)

		T = ThrustTotal(outputs, sections, rotor, op)
		Q = TorqueTotal(outputs, sections, rotor, op)

		return [T; Q]
	end

	#ForwardDiff.jacobian(BEMTwrapper, x)
end 

function eval_f(x)
	n = Nparam
	return - BEMTwrapper(x)[2] * x[3*n+4]
end

function eval_g(x, g)
	T,Q = BEMTwrapper(x)
	g[1] = Q
	return g[2] = T
end

function eval_grad_f(x, grad_f)
    ForwardDiff.gradient!(grad_f, eval_f, x)
end

function eval_jac_g(
    x, 
    rows,
    cols,
    values
)

    if values === nothing
        # m = num of constraint
        # n = num of parameters
        id = 1
        for m in 1:2
            for n in 1:3*Nparam +5
                rows[id] = n
                cols[id] = m
                id += 1
            end
        end
    else
        jac = zeros(2,3*Nparam +5)
        ForwardDiff.jacobian!(jac, (y,x)->eval_g(x,y), zeros(2), x)
        values .= vec( jac ) 
    end
    return nothing
end


#hes = ForwardDiff.jacobian(x -> ForwardDiff.jacobian((x)->begin
#	T,Q = BEMTwrapper(x)
#	[Q; T]
#	end, x), x)

nzJ = 2 *(Nparam*3 +5) 
nzH = 2 *(Nparam*3 +5)^2
n = 3*Nparam +5
m = 2
#x = [r; chord; theta; Rhub; Rtip; Vinf; Omega; rho]
begin

	rP = LinRange( (1 -Rhub/Rtip)/(Ncut+1) +Rhub/Rtip, 1 -(1 -Rhub/Rtip)/(Ncut+1), Nparam ) *Rtip
	chordP_L = 0.03 *ones(Nparam)
	chordP_U = 0.2 *ones(Nparam) *Rtip
	thetaP_L = -pi/6 *ones(Nparam)
	thetaP_U = 5 *ones(Nparam)

	Omega_L = 0.01
	Omega_U = 1220.

	x_L = [rP; chordP_L; thetaP_L; Rhub; Rtip; Vinf; Omega_L; rho]
	x_U = [rP; chordP_U; thetaP_U; Rhub; Rtip; Vinf; Omega_U; rho]
	
end
#g_L = [3.27-0.11; -Inf]
#g_U = [3.27+0.11; 35]

g_L = [3.; 0.01]
g_U = [15.; 25]

prob = Ipopt.CreateIpoptProblem(
    n,
    x_L,
    x_U,
    m,
    g_L,
    g_U,
    nzJ,
    nzH,
    eval_f,
    eval_g,
    eval_grad_f,
    eval_jac_g,
    nothing
)

AddIpoptStrOption(prob, "hessian_approximation", "limited-memory")

prob.x = x

IpoptSolve(prob)

res = prob.x

objective = prob.obj_val

BEMTwrapper(res)

println(res)
println(objective)

"""
    writedlm("BEMTOptiTestRUN.dat", res)

    outputs = BEMTwrapper(res)

    println("BEMT")
    println(BEMTwrapper(res))

    using CCBlade
    function CCbladewrapper(x)
        n = Ncut
        r = x[1:n]
        chord = x[n+1:2*n]
        theta = - x[2*n+1:3*n]

        Rhub = x[3*n+1]
        Rtip = x[3*n+2]

        Vinf = x[3*n+3]
        Omega = x[3*n+4]
        rho = x[3*n+5]

        af = CCBlade.AlphaAF("SG6043_360_Polar_NREL_Format.dat")
        
        rotor = CCBlade.Rotor(Rhub, Rtip, 3; turbine=true, tip=nothing)
        sections = CCBlade.Section.(r, chord, theta, Ref(af))
        op = CCBlade.simple_op.(Vinf, Omega, r, rho)

        sols = CCBlade.solve.(Ref(rotor), sections, op)

        T,Q = thrusttorque(rotor, sections, sols)
        return [T; Q]
    end


    println("CCBlade")
    println(CCbladewrapper(res))



    #x = readdlm("BEMTOptiTestRUN.dat")
    function BEMTwrapper2(x; Ncut=20)
        n = Ncut
        r = x[1:n]
        chord = x[n+1:2*n]
        theta = x[2*n+1:3*n]

        Rhub = x[3*n+1]
        Rtip = x[3*n+2]

        Vinf = x[3*n+3]
        Omega = x[3*n+4]
        rho = x[3*n+5]

        r = LinRange( Rhub/Rtip, 1., Ncut+2 )[2:end-1] *Rtip
        chord = akima(rP, chordP, r)
        theta = akima(rP, thetaP, r)

        rotor = Rotor(Rtip, Rhub, N)
        sections = Section.(r, chord, theta, interpCLCDstruct(airfoildata)..., (r)->1, N*chord ./(2*pi*r))
        op = OpCond(Vinf, Omega, rho)

        sols = Solve.(Ref(rotor), sections, Ref(op))
        outputs = StructArray(sols)

        T = ThrustTotal(outputs, sections, rotor, op)
        Q = TorqueTotal(outputs, sections, rotor, op)

        return [T; Q], outputs
    end


    sols, outs = BEMTwrapper2(OptiRes)

    #outs.Qprime |> plot
"""