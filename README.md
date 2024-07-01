# ItuRPropagations
[![Build Status](https://github.com/disberd/ItuRPropagation.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/disberd/ItuRPropagation.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/disberd/ItuRPropagation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/disberd/ItuRPropagation.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
A Julia implementation of the ITU-Recommendations for space links covering cloud, gaseous, rain, and scintillation attenuations.

> [!NOTE]
> This is a fork of the original repository aimed at adding some functionality and reducing allocations to speed up the computations. This version will be for the time being increased in version and tracked on a private registry for internal use. Inclusion of the changes in this repo to the upstream repository will be done if the original author @HillaryKChao agrees.

> [!WARNING]
> This fork also includes a force push to delete the big files used as input and previously stored inside the `src/data` folder from the git history. This to avoid downloading >100Mb of data at each git clone of the package

## Installation
This fork is not currently registered in the general registry (while the original repository is).
To add it, you have then to explicitly point to this repository with the folloing command in the `Pkg` repl mode
```
add https://github.com/disberd/ItuRPropagation.jl
```
You can check if installation was successful by exiting the package manager and using
```
using ItuRPropagation
```

## ITU-R Recommendations
The following ITU-R Recommendations are implemented:
### Cloud Attenuation
*   **ITU-R P.840-8:** Attenuation due to clouds and fog 
### Gaseous Attenuation
*   **ITU-R P.676-13:** Attenuation by atmospheric gases
*   **ITU-R P.453-14:** The radio refractive index: its formula and refractivity data
*   **ITU-R P.835-16:** Reference Standard Atmospheres
*   **ITU-R P.1511-2:** Topography for Earth-to-space propagation modelling
*   **ITU-R P.2145-0:** Digital maps related to the calculation of gaseous attenuation and related effects
### Rain Attenuation
*   **ITU-R P.618-13:** Propagation data and prediction methods required for the design of Earth-space telecommunication systems
*   **ITU-R P.837-7:** Characteristics of precipitation for propagation modelling
*   **ITU-R P.838-8:** Specific attenuation model for rain for use in prediction methods
*   **ITU-R P.839-4:** Rain height model for prediction methods.
*   **ITU-R P.1511-2:** Topography for Earth-to-space propagation modelling
### Scintillation Attenuation
*   **ITU-R P.618-13:** Propagation data and prediction methods required for the design of Earth-space telecommunication systems
*   **ITU-R P.453-14:** The radio refractive index: its formula and refractivity data

### Noise temperature
*   **ITU-R P.372-16:** Radio noise (only for earth stations - for space stations, a different estimate was used)
*   **ITU-R P.2145-0:** Digital maps related to the calculation of gaseous attenuation and related effects

##  Validation
This implementation has been validated using the [ITU Validation examples (rev 7.0.5)](https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev7.0.5.xlsx).
