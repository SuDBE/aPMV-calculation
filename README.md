# aPMV-calculation
An improved method for calculating aPMV index.

Folder Input includes a cleaned dataset based on the Boxplot rule (python code: remove NA - outliter.py). Data were extracted from the ASHRAE Global Thermal Comfort Database II: Parkinson, Thomas et al. (2022). ASHRAE global database of thermal comfort field measurements [Dataset]. Dryad. https://doi.org/10.6078/D1F671

The PMV calculation can be done by using the package pythermalcomfort, which is developed by: Tartarini, F., Schiavon, S., 2020. pythermalcomfort: A Python package for thermal comfort research. SoftwareX 12, 100578. https://doi.org/10.1016/j.softx.2020.100578

## For ASHRAE data
lambda_N_ASHRAE.py will calculate the λ value of climate Csb from the ASHRAE database as an example.

## For general usage  
lambda_N.py provides a general calculation process for the csv file, which should meet the following requirements:
### Packages 
pandas, numpy, sklearn, matplotlib, seaborn, pythermalcomfort  
### File name  
aPMV.csv  
### File columns (no unit)  
TSV: thermal sensation vote, -3 to +3  
Clo: clothing level, clo  
Met: metabolic rate, met  
Ta: air temperature, ℃  
Tr: radiant temperature, ℃  
RH: relative humidity, %  
Vel: air velocity, m/s    
#### Example:  
![image](https://github.com/SuDBE/aPMV-calculation/assets/68046449/4b0f2998-f5b1-41f5-8f88-70ff26d6b51a)  

### Outputs
#### 1. aPMV curve with λ values
![pic-aPMV curve](https://github.com/SuDBE/aPMV-calculation/assets/68046449/80841f34-5d6a-40ed-8830-6d305f3d4dd8)

#### 2. aPMV predictions compared with TSV and PMV  
![pic-aPMV points](https://github.com/SuDBE/aPMV-calculation/assets/68046449/e95d4a39-b0d9-4080-9305-f5e181c314f8)
