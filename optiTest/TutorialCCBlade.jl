import Pkg; Pkg.activate(".")
using CCBlade

Rtip =  # inches to meters
Rhub = 0.10*Rtip
B = 2  # number of blades

rotor = Rotor(Rhub, Rtip, B)#; turbine=true)#, tip=nothing)

###############

propgeom = [
    0.135/0.55   0.20   -30
    1.00   0.15   -9
    ]

r = LinRange( propgeom[1,1], propgeom[2,1], 100 ) *Rtip
chord = LinRange( propgeom[1,2], propgeom[2,2], 100 ) *Rtip
theta = LinRange( propgeom[1,3], propgeom[2,3], 100 )*pi/180

###############

af = AlphaAF("naca4412.dat")

sections = Section.(r, chord, theta, af)

#########################

Vinf = 10.0
Omega = 12/60*2*pi  # convert to rad/s
rho = 1.225

op = simple_op.(Vinf, Omega, r, rho)

out = solve.(Ref(rotor), sections, op)

out.alpha*180/pi
out.phi
