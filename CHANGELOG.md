# Changelog

This file contains the changelog for the ItuRPropagation package. It follows the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.

## Unreleased

### Changed
- The code for the `IturP835` module has been updated to follow version 7 of the ITU-R P.835 recommendation.
  - The `standardwatervapordensity` function now only accepts the geometric height `Z` as sole positional argument. It is still accepting specific values of `ρ₀` (also as `rho_0`), `T`, and `P` as keyword arguments but they all have default values as per equations 6-8 in Section 1.2 of the recommendation.
