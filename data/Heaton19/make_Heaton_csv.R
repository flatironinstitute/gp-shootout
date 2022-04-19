# R code to convert Heaton '19 comparison datasets to CSV format

load("AllSatelliteTemps.RData")
load("SatelliteTemps.RData")
load("AllSimulatedTemps.RData")
load("SimulatedTemps.RData")

# At this point you get:
# > ls()
# [1] "all.sat.temps" "all.sim.data"  "sat.temps"     "sim.data"

write.csv(all.sat.temps,"AllSatelliteTemps.csv")
write.csv(sat.temps,"SatelliteTemps.csv")
write.csv(all.sim.data,"AllSimulatedTemps.csv")
write.csv(sim.data,"SimulatedTemps.csv")
