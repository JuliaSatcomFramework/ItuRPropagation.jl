# Changelog

This file contains the changelog for the ItuRPropagation package. It follows the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.

## Unreleased

### Changed
- All instances of the word `vapor` in the codebase has been replaced with `vapour` for consistency to the Recommendation texts
- The code for the `IturP835` module has been updated to follow version 7 of the ITU-R P.835 recommendation.
  - The `standardwatervapourdensity` function now only accepts the geometric height `Z` as sole positional argument. It is still accepting specific values of `ρ₀` (also as `rho_0`), `T`, and `P` as keyword arguments but they all have default values as per equations 6-8 in Section 1.2 of the recommendation.
- The `ItuRP2145` module has been refactored to use an independent artifact and for a more modular structure.

### Added
- The `ItuRP2145` now supports computation of average annual values, as opposed to only the CCDF values supported in previous versions
  - The module now also supports computing the total barometric pressure and the integrated water vapour content, available via the `surfacepressureannual` and `surfacewatervapourcontentannual` functions, respectively.
- The function of the `ItuRP2145` module are now tested against the validation examples excel provided by ITU-R.
