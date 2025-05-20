# Nonergodic_adjustments_CEUS
###########################################################################################
# Nonergodic site, source and path adjustments for Central and Eastern United States (CEUS)
###########################################################################################

This repository contains the results of nonergodic site, source, and path adjustments for the Central and Eastern United States (CEUS), developed by Mehran Davatgari-Tafreshi and Shahram Pezeshk (2025).

######################
## Folder Structure
######################

Source_Mean.csv: This file contains the nonergodic source adjustments for PGA and for periods ranging from 0.01 to 10 seconds (22 points). Latitude and longitude in this file indicate the coordinates of each source region (at the midpoint).

Source_Epistemic.csv: This file contains the epistemic uncertainty associated with the nonergodic source adjustments (standard error of the mean) for PGA and for periods ranging from 0.01 to 10 seconds (22 points). Latitude and longitude in this file indicate the coordinates of each source region (at the midpoint).

Mean_Site.csv:This file contains the nonergodic site adjustments for PGA and for periods ranging from 0.01 to 10 seconds (22 points).StationLatitude and StationLongitude in this file indicate the coordinates of each station.
StationID, Vs30, and CoastalPlainsSedThick refer to the station name, Vs30 (the time-averaged shear-wave velocity within the top 30 meters of the Earth's surface, in m/s), and the sediment thickness (sediment deposits beneath the Atlantic and Gulf Coastal Plains, in meters), respectively.

Epistemic_Site.csv:This file contains the epistemic uncertainty associated with the nonergodic site adjustments (standard error of the mean) for PGA and for periods ranging from 0.01 to 10 seconds (22 points). StationLatitude and StationLongitude in this file indicate the coordinates of each station. StationID, Vs30, and CoastalPlainsSedThick refer to the station name, Vs30 (the time-averaged shear-wave velocity within the top 30 meters of the Earth's surface, in m/s), and the sediment thickness (sediment deposits beneath the Atlantic and Gulf Coastal Plains, in meters), respectively.

Path_Mean.csv:This file contains the nonergodic path adjustments for each source–site pair for PGA and for periods ranging from 0.01 to 10 seconds (22 points). StationLatitude and StationLongitude indicate the coordinates of each station, while EarthquakeLatitude and EarthquakeLongitude indicate the coordinates of each earthquake.PreferredMagnitude, EarthquakeId, RuptureDistance, StationID, Vs30, and CoastalPlainsSedThick refer, respectively, to the moment magnitude of each earthquake, earthquake ID, rupture distance, station name, Vs30 (the time-averaged shear-wave velocity within the top 30 meters of the Earth's surface, in m/s), and sediment thickness (sediment deposits beneath the Atlantic and Gulf Coastal Plains, in meters).

Path_Epistemic2.csv:This file contains the epistemic uncertainty associated with the nonergodic path adjustments (standard error of the mean) for each source–site pair for PGA and for periods ranging from 0.01 to 10 seconds (22 points). StationLatitude and StationLongitude indicate the coordinates of each station, while EarthquakeLatitude and EarthquakeLongitude indicate the coordinates of each earthquake.PreferredMagnitude, EarthquakeId, RuptureDistance, StationID, Vs30, and CoastalPlainsSedThick refer, respectively, to the moment magnitude of each earthquake, earthquake ID, rupture distance, station name, Vs30 (the time-averaged shear-wave velocity within the top 30 meters of the Earth's surface, in m/s), and sediment thickness (sediment deposits beneath the Atlantic and Gulf Coastal Plains, in meters).

Coefficients.txt: This file contains the coefficients (c1–c11) of the ground motion model developed by Pezeshk et al. (2018) for areas outside the Coastal Plain (see Table 4 of Pezeshk et al., 2018). Each column presents the values of c1–c11 for PGA and for periods ranging from 0.01 to 10 seconds (22 points).

vector.csv:This file contains the coefficients (c1–c5) of the ground motion model developed by Akhani et al. (2024) for areas inside the Coastal Plain (see Table 1 of Akhani et al. 2024). Each column presents the values of c1–c5 for PGA and for periods ranging from 0.01 to 10 seconds (22 points).

lon_lat_cp.txt:This file contains the latitude and longitude coordinates defining the boundary of the Coastal Plain (Boyd et al. 2024).

################
## References
################

Boyd, O. S., Churchwell, D., Moschetti, M. P., Thompson, E. M., Chapman, M. C., Ilhan, O., ... and Rezaeian, S. (2024). Sediment thickness map of United States Atlantic and Gulf Coastal Plain Strata, and their influence on earthquake ground motions. Earthquake Spectra, 40(1), 89-112. DOI: https: 10.1177/87552930231204880.

Pezeshk, S., Zandieh, A., Campbell, K. W., and Tavakoli, B. (2018). Ground‐motion prediction equations for central and eastern North America using the hybrid empirical method and NGA‐West2 empirical ground‐motion models. Bulletin of the Seismological Society of America, 108(4), 2278–2304. DOI: 10.1785/0120170179.
 
Akhani, M., Davatgari-Tafreshi, M., and Pezeshk, S. (2024). Adjusting Central and Eastern United States ground-motion models for use in the Coastal Plain considering the sediment thickness. Earthquake Spectra, 40(4), 2669–2691. DOI: 10.1177/87552930241258354.
