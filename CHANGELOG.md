# Changelog

This file contains the changelog for the ItuRPropagation package. It follows the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.

## Unreleased

Note that basically all changes above are **BREAKING**

### Changed
- All instances of the word `vapor` in the codebase has been replaced with `vapour` for consistency to the Recommendation texts.
- The code for the `IturP835` module has been updated to follow version 7 of the ITU-R P.835 recommendation.
  - The `standardwatervapourdensity` function now only accepts the geometric height `Z` as sole positional argument. It is still accepting specific values of `ρ₀` (also as `rho_0`), `T`, and `P` as keyword arguments but they all have default values as per equations 6-8 in Section 1.2 of the recommendation.
- The `ItuRP2145` module has been refactored to use an independent artifact and for a more modular structure.
- The implementation of module `ItuRP1511` has been refactored to use bicubic interpolation from the `ItuRP1144` module and to rely on different indpendent artifact for 1511 data. It also now updated to version 3 of ITU-R P.1511
- Refactored the internal implementation of the `ItuRP676` to rely on some artifacts for Part1 and Part2 annexes rather than hardcoding in file.
  - The tests have also been changed to directly parse the data from the ITU validation excel
- The `ItuRP676.gaseousattenuation` is now correctly respecting the implementation in Annex 2 of ITU-R P.676-13, computing statistical gaseous attenuation between 1 and 350GHz
- The implementation of cloud attenuation in module `ItuRP840` has been refactored, now relying on artifacts and on convenience interpolation structures in `ItuRP1144`
- The `ItuRP840.cloudattenuation` now has a different signature, with the elevation and exceedance probability swapped to e more consistent with the recent changes to other module.
- Refactored the `ItuRP453` module to use independent artifacts and now provide a function to interpolate the wetterm refractive index not only at the 50 percentile
- Refactored the `ItuRP839` module to use specific artifacts and its test to use the ITU validation excel version 8.3

### Added
- A new `ItuRP1144` module has been added to hold the interpolation functions defined recommendation ITU-R P.1144-12.
  - For the moment, this only implements the bilinear interpolation on a square grid and bicubic interpolation.
  - Added structurs for convenience interpolation with bilinear square grid
- The `ItuRP2145` module now supports computation of average annual values, as opposed to only the CCDF values supported in previous versions
  - The module now also supports computing the total barometric pressure and the integrated water vapour content, available via the `surfacepressureannual` and `surfacewatervapourcontentannual` functions, respectively.
- The function of the `ItuRP2145` module are now tested against the validation examples excel provided by ITU-R.
- Added a constant `SUPPRESS_WARNINGS::Ref{Bool}` in the main module that can be used to suppress warnings given by the various functions (Defaults to `false`)