import Pkg; Pkg.activate(".")
using CCBlade

Rtip = 0.55
Rhub = 0.135
N = 3
Vinf = 10.
Ncut = 100
Omega = 500/60*2*pi
rho = 1.225
propgeom = [0.135/0.55   0.20   -30
            1.00         0.15   -9]

r = LinRange( Rhub/Rtip, 1., Ncut+2 )[2:end-1] *Rtip
chord = LinRange( propgeom[1,2], propgeom[2,2], Ncut+2 )[2:end-1] *Rtip;           
theta = LinRange( propgeom[1,3], propgeom[2,3], Ncut+2 )[2:end-1] *pi/180;  

rotor = CCBlade.Rotor(Rhub, Rtip, B; turbine=true, tip=nothing)


af = AlphaAF("SG6043_360_Polar_NREL_Format.dat")

sections = CCBlade.Section.(r, chord, theta, Ref(af))


Vinf = 10.0
Omega = 500/60*2*pi  # convert to rad/s
rho = 1.225

op = simple_op.(Vinf, Omega, r, rho)

out = CCBlade.solve.(Ref(rotor), sections, op)

thrusttorque(rotor, sections, out)


out.alpha*180/pi
out.phi

