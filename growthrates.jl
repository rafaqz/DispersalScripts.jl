using GrowthRates, Dates
using Unitful: °C, K, hr, d
using CuArrays

import GrowthRates: arraysetup
arraysetup(a) = CuArray(a)
# arraysetup(a) = a

#################################################################
# User Script

p25 = 3.377850e-01
H_A = 3.574560e+04K
H_L = -1.108990e+05K
T_0_5L = 2.359187e+02K
H_H = 3.276604e+05K
T_0_5H = 2.991132e+02K
R = 1.987

ipgrowth = IntrinsicPopGrowth(p25, H_A, H_L, T_0_5L, H_H, T_0_5H, 1/K)

lowersm = 0.0 # no data
uppersm = 1.0 # no data
lowerct = 7.0°C |> K  # Enriquez2017
lowerctm = -log(1.00) * K^-1
upperct = 30.0°C |> K # Kimura2004
upperctm = -log(1.15) * K^-1
lowerwilt = 0.5 # default?
wiltingm = -log(1.1)

# general options
cold = ColdDays(lowerct, lowerctm)
hot = HotDays(upperct, upperctm)
wilt = Wilting(lowerwilt, wiltingm)

start_date = Date("2016-01-01")
end_date = Date("2016-12-31")
# end_date = Date("2016-02-28")
data_path = "/home/raf/Storage/drosophila/climate/SMAP_raw/" # raw climatic data
source = SMAP()
models = (ipgrowth, cold, hot, wilt)

stresses = @time growthrates(models, data_path, start_date, end_date, 365.25d/12)
