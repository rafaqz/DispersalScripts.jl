using HDF5, GrowthRates, GeoData, Dates, Plots
using Unitful: °C, K, hr, d
using GeoData: Time

series = SMAPseries("smap"; lazyview=(Lon<|Between(-125, -75), Lat(150:480)))
series = series[Between(Date(2016), Date(2018))]

month = Second(356.25*24*60*60/12) 
period = month
sampleperiod = Day(1)

t = DateTime(2012) + month

round(t, Day(1))


# parameters
p25 = 3.377850e-01
H_A = 3.574560e+04K
H_L = -1.108990e+05K
T_0_5L = 2.359187e+02K
H_H = 3.276604e+05K
T_0_5H = 2.991132e+02K
R = 1.987

lowersm = 0.0 # no data
uppersm = 1.0 # no data
lowerct = 7.0°C |> K  # Enriquez2017
lowerctm = -log(1.00) * K^-1
upperct = 30.0°C |> K # Kimura2004
upperctm = -log(1.15) * K^-1
lowerwilt = 0.5 # default?
wiltingm = -log(1.1)

# models
growth = IntrinsicPopGrowth{:surface_temp}(1/K, p25, H_A, H_L, T_0_5L, H_H, T_0_5H)
cold = ColdDays{:surface_temp}(lowerct, lowerctm)
hot = HotDays{:surface_temp}(upperct, upperctm)
wilt = Wilting{:land_fraction_wilting}(lowerwilt, wiltingm)
models = (growth, cold, hot, wilt)
# models = (growth, wilt)

@time output = growthrates(models, series, timestep, sumsteps)
output[1] |> plot
