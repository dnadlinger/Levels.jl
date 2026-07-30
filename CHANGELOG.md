# CHANGELOG

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog],
and this project adheres to [Semantic Versioning].

## [Unreleased]

- Initial release
- `StateBasis` for ordered Zeeman-state index bookkeeping, plus the generic
  `state_pairs` transition enumeration between two levels
- Zeeman structure: `lande_g` (LS coupling with reduced-mass correction, using
  the new species `mass` field), `zeeman_shift`, `zeeman_sensitivity`, and
  `zeeman_hamiltonian` for arbitrary field directions
- Electric-quadrupole coupling geometry: `beam_vectors`, `quadrupole_geometry`,
  `quadrupole_couplings`, and `rabi_normalised`
- ac Stark (light) shifts of states and transitions: `light_shift`, with the
  atomic structure hoisted into a reusable `LightShiftCoefficients`, plus the
  `scalar_polarisability`/`vector_polarisability`/`tensor_polarisability`
  decomposition; driven by the new species `polarisabilities` field, populated
  for the ⁸⁸Sr⁺ S``_{1/2}`` and D``_{5/2}`` clock levels
- Near-resonant electric-quadrupole light shifts — the shift a probed
  quadrupole Zeeman component receives from the components sharing one of its
  states, the "E2 ac Stark shift" of optical-clock evaluations:
  `quadrupole_light_shift`, also folded into `light_shift` via the `n`/`B`
  keyword arguments, with the per-channel geometry weights exposed as
  `Levels.quadrupole_weights`
- `Levels.PeriodicDynamics` submodule for transitions under periodic driving
  (trap rf): the `DrivenTransition` rotating-frame model with ac Zeeman drives
  and laser phase modulation, a dressed-state Floquet engine
  (`dress`/`sideband_amplitude`), and an exact monodromy-matrix engine
  (`exact_sideband`/`stroboscopic_populations`) as nonperturbative cross-check

<!-- Links -->

[keep a changelog]: https://keepachangelog.com/en/1.1.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

<!-- Versions -->

[unreleased]: https://github.com/dnadlinger/Levels.jl/compare/v0.1.0...HEAD
