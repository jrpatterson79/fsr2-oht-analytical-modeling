## Overview
The data and code in this repository is used to generate the analysis presented in Patterson and Cardiff (2023), _Do Simple Analytical Models Capture Complex Fracture Hydraulics? Oscillatory Flow Tests Suggest Not_ This paper describes an oscillatory hydraulic tomography campaign conducted in 2019 near Madison, WI and uses simple analytical modeling approaches to characterize the hydraulic properites of a single bedrock fracture and attempts to explore the mechanisms leading to period-dependent parameters under oscillatory flow conditions. 

## Field Data
This folder contains all of the raw field data collected during oscillatory flow experiments conducted at (FSR)^2 in 2019. It also contains all of the output .mat files generated during data processing and quality control. Though all of the data processing code is provided with the supplementary information, it is not necessary to go through these steps to repeat the analysis generated in this paper.
* 2019 Test Database.xlsx contains all of the oscillatory flow tests conducted in 2019 with the relevant testing parameters and packer depths relative to top of casing. Rows highlighted with the same color indicate repeat tests. This file is needed to run some of the data processing / analysis code described below.
* FSR2_WellSurvey.csv contains the geodetic survey of the (FSR)^2 wells conducted in 2019. Data was collected using EMLID RTK units. This file is needed to run some of the data processing / analysis code described below.

## Analysis Code
At a high level, the code in this repository takes raw signals as measured during oscillatory flow testing and implements a data processing workflow - with intermediate quality control steps - to estimate the hydraulic properties - with uncertainty - of a single bedrock fracture using a gradient-based inversion strategy described in Patterson and Cardiff (2022). The analysis presented in this paper was generated using code developed and executed in MATLAB 2019b. At the time of upload the code has been tested and runs without issues in MATLAB 2021b. 

* validation_2D.m validates the numerical model used for the analysis in the **Borehole Storage** section.
* borehole-storage.m is a 2D numerical model used to conduct the analysis found in the **Borehole Storage** section.
* All .m files found in the Func_Lib subdirectory are dependencies that are called in the main codes mentioned above. Documentation and variable description is found at the top of each individual file.

### Field Data Analysis
* step1.m - step7.m conducts the field data analysis described in the **Field Analysis** and **Parameter Estimation** sections. The necessary .mat files are provided, so that steps 6b, 7, and 7b can all be run without repeating our data processing in steps 1-5.

## License
The code and data are provided as open source under the GNU General Public License v3.0. It is provided without warranty, but should perform as described in the manuscript when executed without modification. If you would like to use any of the code in this repository for research, software, or publications, we ask that you provide a citation to the code and journal article (See references below).

## References
Patterson, Jeremy R., and M. Cardiff. 2022. Aquifer Characterization and Uncertainty in Multi-Frequency Oscillatory Flow Tests: Approach and Insights. Groundwater 60, no. 2: 180–191, https://doi.org/10.1111/gwat.13134.

Patterson, Jeremy R., and M. Cardiff. 2023. Do Simple Analytical Models Capture Complex Frature Hydraulics? Oscillatory Flow Tests Suggest Not. Groundwater, https://doi.org/10.1111/gwat.13297