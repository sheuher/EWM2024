import Pkg; Pkg.activate(".")
using CSV, DelimitedFiles, DataFrames
using Plots
using Distributions, Statistics

windDatCSV = "/home/sheuher/sciebo - Tey, Sheu Her (soshteyy@uni-duisburg-essen.de)@uni-duisburg-essen.sciebo.de/EWM_groupA/PCRoomBackup2024MAI22/02_Randbedingungen_Winddaten/Winddaten/Windgeschwindigkeit - Keschhütte.csv"
isfile(windDatCSV)

windDF = CSV.read(windDatCSV, DataFrame)

function comma2dot(x::AbstractString)
    x = replace(x, ","=>".")
    parse(Float64, x)
end

windSpeeds = comma2dot.(windDF[:,3])

windSpeeds = windSpeeds[windSpeeds .!= 0.]

windSpeedsDist =  fit_mle(Weibull, windSpeeds)

begin

histogram(windSpeeds, bins=0:0.75:40, normalize=:probability)
plot!(0:0.1:40, pdf.(windSpeedsDist, 0:0.1:40), c=:red)

end