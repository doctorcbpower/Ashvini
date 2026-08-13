# Ashvini: physical models and parameter reference

This document lists the governing equations behind every module in [`ashvini/`](ashvini/) and every parameter that controls them, with defaults as shipped in [`run_params.yaml`](run_params.yaml). For module-to-file mapping and usage, see [README.md](README.md). Line-level detail (exact closed-form solutions used by the vectorised fast path) is in the module docstrings/comments; this document focuses on the physics and the config surface.

Units: masses in Msun, time in Gyr, metallicities in Z_solar units, unless noted. Two different `astropy` cosmologies are used internally: [`utils.py`](ashvini/utils.py) (cosmic time/redshift interpolation, used everywhere) uses **Planck18**; [`reionization.py`](ashvini/reionization.py) uses **Planck15** for its own `H(z)` term. This is a pre-existing inconsistency between modules, not a deliberate modelling choice -- worth knowing if you're chasing small discrepancies in reionization suppression near the LCDM/`H(z)` boundary.

Every reservoir (gas, stars, gas metals, stellar metals, dust, black hole) is evolved as an ODE in cosmic time along a halo's merger history. `main.py` provides two integrators over the same equations: `run1_scalar()` (per-halo, per-timestep `scipy.integrate.solve_ivp`, kept as a trusted reference) and `run_forest()`/`run1()` (the production path, vectorised across all haloes in a forest, using closed-form/quadrature updates instead of `solve_ivp` -- see [Numerical scheme](#numerical-scheme) at the end). Both solve *exactly* the same equations below; they differ only in how the ODEs are stepped.

## Free-fall time

Used by both star formation and BH growth:

```
t_ff(z) = 0.141 * t_Hubble(z)     [Gyr],   t_Hubble(z) = 1/H(z)
```

## Star formation ([`star_formation.py`](ashvini/star_formation.py))

```
dM_star/dt = SFR(t, M_gas) = (epsilon_ff / t_ff(z(t))) * M_gas
```

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| epsilon_ff | `star_formation.efficiency` | 0.015 | Star formation efficiency per free-fall time |

## Gas reservoir ([`gas_evolve.py`](ashvini/gas_evolve.py))

Cosmological inflow, suppressed by the UV background if enabled:

```
Mdot_gas,inflow(z, M_halo, Mdot_halo) = (Omega_b/Omega_m) * Mdot_halo * S_uv(z, M_halo, Mdot_halo)
```

Gas mass ODE (all baryonic sinks acting on the reservoir):

```
dM_gas/dt = Mdot_acc - SFR(t,M_gas) - eta_w(z,M_halo,Z_star)*SFR_wind - Mdot_BH - Mdot_AGN,wind
```

- `eta_w` is the supernova mass-loading factor ([below](#supernova-feedback-supernovae_feedbackpy)).
- `SFR_wind` is the SFR driving the SN wind, selected by `supernova.type`:
  - `"no"` -> `SFR_wind = 0` (feedback off)
  - `"instantaneous"` -> `SFR_wind = SFR(t, M_gas)` (present-time SFR)
  - `"delayed"` (default) -> `SFR_wind` = the SFR from `supernova.delay_time` Gyr ago
- `Mdot_BH` is instantaneous BH accretion -- gas flowing into the BH leaves the reservoir the moment it's accreted.
- `Mdot_AGN,wind = eta_agn * Mdot_BH_delayed` ([AGN feedback](#agn-feedback-agn_feedbackpy)).

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| UV_background enabled | `reionization.UVB_enabled` | True | Whether inflow is suppressed by the UV background |

## Metal enrichment ([`metallicity.py`](ashvini/metallicity.py))

Gas-phase metal mass, with a separate branch for halos with no gas (pure inflow + delayed-wind enrichment, no dilution term):

```
M_gas = 0:   dM_Z,gas/dt = Z_IGM*Mdot_acc + Z_yield*SFR_wind

M_gas > 0:   dM_Z,gas/dt = Z_IGM*Mdot_acc - M_Z,gas*SFR(t,M_gas)/M_gas + Z_yield*SFR_wind
                            - eta_w(z,M_halo,Z_star) * (M_Z,gas/M_gas) * SFR_wind
```

Stellar metal mass (new stars inherit the gas-phase metal fraction at formation):

```
dM_Z,star/dt = SFR(t,M_gas) * M_Z,gas/M_gas     (0 if M_gas = 0)
```

Stellar metallicity `Z_star = M_Z,star/M_star` (0 if `M_star=0`) feeds back into the mass-loading factor and dust wind-loss below.

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| Z_IGM | `metallicity.Z_IGM` | 1.0e-3 | Metallicity of inflowing IGM gas (Z_solar) |
| Z_yield | `metallicity.Z_yield` | 0.06 | Metal yield per unit (delayed) SFR returned to the ISM |

## Dust ([`dust.py`](ashvini/dust.py))

```
growth_rate      = Y_d * SFR_past
SNe_rate         = Gamma * SFR_past / M_star,past          (0 if M_star,past <= 0)
dust_loading      = 1 - exp(-(M_gas/M_crit)^Alpha)           (saturates to 1 for M_gas >> M_crit)
destruction_rate  = M_swept * SNe_rate * dust_loading
wind_loss         = (M_dust/M_gas) * eta_w(z,M_halo,Z_star) * SFR_past    (0 if M_gas <= 0)

dM_dust/dt = growth_rate - destruction_rate - wind_loss
```

`SFR_past` and `M_star,past` are the SFR and the *stellar mass formed in that single delayed step interval* (not cumulative) from `supernova.delay_time` ago -- this represents the population of young SN progenitors formed at that specific past episode, driving dust destruction now. Both are 0 until enough cosmic time has elapsed to have delay history.

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| M_swept | `dust.m_swept` | 1.0e3 Msun | ISM mass swept/destroyed per SN (`= epsilon_d * M_ISM`) |
| Y_d | `dust.dust_yield` | 0.004 | Dust mass produced per unit (delayed) SFR |
| Gamma | `dust.dust_gamma` | 1.3e-4 | SN rate per unit stellar mass formed in the delayed episode |
| Alpha | `dust.dust_alpha` | 8 | Steepness of the destruction-efficiency saturation curve |
| M_crit | `dust.m_crit` | 1.0e5 Msun | Gas-mass scale at which destruction efficiency saturates |

## Supernova feedback ([`supernovae_feedback.py`](ashvini/supernovae_feedback.py))

Metallicity modulation (smooth logistic step from strong feedback at low `Z_star` to weak feedback at high `Z_star`; the step center/width/floor/ceiling below are hard-coded, not exposed in the YAML):

```
f_Z(Z_star) = sigmoid(-(Z_star - 0.1)/0.01) * (0.25 - 1) + 1
```

Mass-loading factor (mass of gas ejected per unit SFR):

```
eta_w(z, M_halo, Z_star) = epsilon_p * pi_fid * (10^11.5 / M_halo)^(1/3) * (9/(1+z))^(1/2) * f_Z(Z_star)
```

`M_halo = 0` (unformed progenitor) returns `eta_w = 0` exactly, avoiding a `1/M_halo` blow-up.

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| type | `supernova.type` | `"delayed"` | `"no"` / `"instantaneous"` / `"delayed"` -- which SFR drives the wind (see [gas reservoir](#gas-reservoir-gas_evolvepy)) |
| delay_time | `supernova.delay_time` | 0.015 Gyr | Lookback time for the delayed-feedback SFR |
| epsilon_p | `supernova.epsilon_p` | 5 | Mass-loading normalization |
| pi_fid | `supernova.pi_fid` | 1 | Fiducial mass-loading efficiency multiplier |

## Black hole growth ([`black_holes_growth.py`](ashvini/black_holes_growth.py))

Eddington rate (mass-proportional):

```
Mdot_Edd(M_BH) = kappa * M_BH,   kappa = 4*pi*G*m_p / (sigma_thomson * c)   [consistent Msun/Gyr units]
```

Growth rate, gas-supply-limited but capped at a (possibly super-Eddington) multiple of the Eddington rate:

```
dM_BH/dt = min( (epsilon_BH / t_ff(z)) * M_gas,  f_Edd * kappa * M_BH )
```

**Seeding** -- three independently enable-able channels, checked in priority order (first satisfied wins per halo, using the *previous* step's halo mass/redshift/gas metallicity), only for not-yet-seeded halos:

| Channel | Config block | Condition | Seed mass |
|---|---|---|---|
| Pop III remnant | `black_holes.seeding.pop3` | `z >= z_min` AND `M_halo >= M_halo_min` AND `Z_gas <= Z_gas_max` | `M_seed` |
| Direct collapse | `black_holes.seeding.direct_collapse` | (as above, checked next) | `M_seed` |
| Halo-mass threshold | `black_holes.seeding.halo_mass_threshold` | `M_halo >= M_halo_min` (fallback) | `M_seed` |

Default values:

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| enabled | `black_holes.seeding.pop3.enabled` | True | Pop III channel on/off |
| z_min | `black_holes.seeding.pop3.z_min` | 15 | Only forms this early in the universe |
| M_halo_min | `black_holes.seeding.pop3.M_halo_min` | 1.0e6 Msun | Minihalo mass threshold (H2-cooling, Tvir ~ 1e4 K) |
| Z_gas_max | `black_holes.seeding.pop3.Z_gas_max` | 1.0e-4 | Near-primordial gas required |
| M_seed | `black_holes.seeding.pop3.M_seed` | 1.0e2 Msun | Light seed (Pop III remnant) |
| enabled | `black_holes.seeding.direct_collapse.enabled` | True | Direct-collapse channel on/off |
| z_min | `black_holes.seeding.direct_collapse.z_min` | 10 | Redshift floor |
| M_halo_min | `black_holes.seeding.direct_collapse.M_halo_min` | 1.0e7 Msun | Atomic-cooling halo threshold |
| Z_gas_max | `black_holes.seeding.direct_collapse.Z_gas_max` | 1.0e-5 | Pristine gas required to suppress H2 fragmentation |
| M_seed | `black_holes.seeding.direct_collapse.M_seed` | 1.0e5 Msun | Heavy seed |
| enabled | `black_holes.seeding.halo_mass_threshold.enabled` | True | Fallback channel on/off |
| M_halo_min | `black_holes.seeding.halo_mass_threshold.M_halo_min` | 1.0e10 Msun | Generic mass threshold |
| M_seed | `black_holes.seeding.halo_mass_threshold.M_seed` | 1.0e3 Msun | Seed mass |
| efficiency | `black_holes.efficiency` | 0.001 | `epsilon_BH`, gas-supply growth efficiency (same functional form as SF efficiency) |
| eddington_multiplier | `black_holes.eddington_multiplier` | 1.0 | `f_Edd`, Eddington cap multiplier (>1 allows super-Eddington growth) |

## AGN feedback ([`agn_feedback.py`](ashvini/agn_feedback.py), [`black_holes_growth.py`](ashvini/black_holes_growth.py))

```
Mdot_AGN,wind = eta_agn_eff * Mdot_BH,delayed
```

`Mdot_BH,delayed` is the BH growth rate from `feedback_delay_time` Gyr ago (0 if no history yet); with the default `feedback_delay_time = 0.0`, this reduces to the current step's own accretion rate (instantaneous feedback). By default `eta_agn_eff` is just the constant `eta_agn`; enabling `black_holes.sigma_feedback` replaces it with the M-sigma self-regulation coupling below.

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| eta_agn | `black_holes.eta_agn` | 0.5 | AGN wind mass-loading efficiency |
| feedback_delay_time | `black_holes.feedback_delay_time` | 0.0 Gyr | Lag between BH accretion and the wind it powers; 0.0 = instantaneous |

### Isothermal-sphere M-sigma self-regulation (optional)

King (2003, 2005) and Power, Zubovas, Nayakshin & King (2011, MNRAS 413, L110) argue AGN wind coupling to bulge gas isn't a constant efficiency: it's weak (momentum-driven, Compton-cooled, effectively trapped near the BH) while `M_BH < M_sigma`, and becomes effective at expelling gas once `M_BH >= M_sigma`, where `M_sigma` is set by balancing the wind's Eddington-limited momentum thrust against the weight of the overlying gas in an isothermal sphere. `black_holes.sigma_feedback.enabled` (default `False`) swaps the constant `eta_agn` above for this mechanism; disabled, behaviour is identical to the constant coupling. Only the isothermal-sphere approximation is implemented -- an NFW version (radius-dependent weight, concentration-mass relation, transcendental `M_sigma` solve) is a deliberate follow-up, not implemented here.

Isothermal-sphere velocity dispersion (`velocity_dispersion`, from the virial relations `R_v = sigma/(5*sqrt(2)*H(z))`, `M_v = 2*sigma^2*R_v/G`, PZNK11 eq. 18):

```
sigma(M_halo, z) = [ (5*sqrt(2)/2) * G * H(z) * M_halo ]^(1/3)
```

`H(z)` uses the same Planck18 cosmology as everywhere else in the package (`utils.Hubble_time`), not a third cosmology alongside the Planck18/Planck15 split noted above. `M_halo <= 0` gives `sigma = 0`.

Self-regulated BH mass scale (`m_sigma`, King 2003/2005, PZNK11 eq. 5):

```
M_sigma(sigma) = (f_g * kappa_es / (pi * G^2)) * sigma^4
```

`kappa_es = sigma_thomson/m_p` is the electron-scattering opacity -- a different quantity from the Eddington-rate constant `kappa` in [Black hole growth](#black-hole-growth-black_holes_growthpy) above, despite the shared name in the literature. `sigma = 0` gives `M_sigma = 0`.

Effective coupling (`coupling_switch`, replacing the constant `eta_agn`):

```
eta_agn_eff(M_BH, M_sigma) = eta_agn * f_switch(M_BH / M_sigma)
f_switch(x) = 1 / (1 + exp(-log10(x) / transition_width))     (logistic in log10(M_BH/M_sigma), centered at M_BH = M_sigma)
```

Neither paper specifies a particular smooth interpolation between the trapped and escaping regimes -- only that the regime change happens near `M_BH ~ M_sigma` -- so `f_switch`'s exact functional form is a deliberate modelling choice, isolated in `coupling_switch()` so it can be revisited independently of the rest of the pipeline. `M_BH <= 0` or `M_sigma <= 0` gives `f_switch = 0` (no coupling without a BH, or for an unformed halo).

**Numerical scheme note**: because `eta_agn_eff` depends on `M_BH(t)`, which evolves within the same step as the gas-mass ODE it enters, `M_BH` is evaluated at the step midpoint (via the same closed-form BH-growth update used for `bh_mass[j]` itself, but over `dt/2`) rather than at the step's own not-yet-known end -- consistent with how every other time-varying coefficient in the gas ODE is already frozen mid-step (see [Numerical scheme](#numerical-scheme)), and necessary to keep the gas-mass ODE affine so `_linear_ode_step` still applies exactly.

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| enabled | `black_holes.sigma_feedback.enabled` | False | Master switch; disabled reproduces the constant `eta_agn` coupling exactly |
| f_g | `black_holes.sigma_feedback.f_g` | 0.16 | Baryon fraction relative to dark matter |
| kappa_es | `black_holes.sigma_feedback.kappa_es` | null (-> `sigma_thomson/m_p`) | Electron-scattering opacity, cm^2/g |
| transition_width | `black_holes.sigma_feedback.transition_width` | 0.1 | Dex width of the smooth `M_BH/M_sigma` switch |

Out of scope for this mechanism (potential follow-ups): an NFW-profile version of `sigma`/`M_sigma`; any feedback-driven cap or quenching of BH growth itself tied to `M_sigma` (PZNK11 eq. 20-23's deviation term); and the competing nuclear-cluster feedback channel (Nayakshin, Wilkinson & King 2009).

**Important caveat, checked directly (see `tests/test_sigma_feedback.py::test_sigma_feedback_wind_alone_does_not_cap_bh_growth_near_m_sigma`): this does not actually hold `M_BH` near `M_sigma`.** "Self-regulation" in King (2003, 2005)/PZNK11 refers to the BH's own growth eventually being throttled once the wind is strong enough to expel/starve the gas supply feeding it -- but this implementation only adds a gas-removal *wind* term (`Mdot_wind = eta_agn_eff * Mdot_BH`); it does not add the explicit growth-quenching term (PZNK11 eq. 20-23) called out as out-of-scope above. Regulation, if any, can only happen indirectly, through the wind depleting `gas_mass` enough to lower the gas-supply-limited growth rate (`epsilon_BH/t_ff * gas_mass`) that feeds `M_BH` itself. Checked over a long (z=25 -> z=0.5), sustained-accretion integration: that indirect loop is nowhere near strong enough here -- cosmological gas resupply vastly outpaces even a maximal wind (`eta_agn=1.0`), so `M_BH` keeps growing over two orders of magnitude past `M_sigma`, not settling near it. In other words, as implemented, `sigma_feedback` is a **wind-strength switch gated by `M_BH/M_sigma`**, not a BH-mass cap -- don't rely on it to keep `M_BH` near the M-sigma relation over cosmic time; that would require the explicitly out-of-scope growth-quenching term above.

## Reionization ([`reionization.py`](ashvini/reionization.py))

Okamoto et al. (2008)-style suppression of baryonic inflow below a characteristic halo mass `M_c(z)`, active only for `z <= 10` (identically 1, i.e. no suppression, above that):

```
M_c(z) = 1.69e10 * exp(-0.63*z) / (1 + exp((z/beta)^gamma))     [Msun]
beta   = z_reion * ( ln(1.82e3 * exp(-0.63*z_reion) - 1) )^(-1/gamma)

S_uv(z, M_halo, Mdot_halo) = max( 0,
    s(M_halo/M_c(z), omega) * [ (1+X) - 2*epsilon(z)*M_halo*X*(1+z)*H(z)/Mdot_halo ] )

s(x,y)  = (1 + (2^(y/3)-1) * x^-y)^(-3/y)
X       = 3*c_omega*(M_c(z)/M_halo)^omega / (1 + c_omega*(M_c(z)/M_halo)^omega),   c_omega = 2^(omega/3)-1
epsilon(z) = d/dz [ ln(1 + exp((z/beta)^gamma)) ]  (evaluated analytically)
```

`M_halo <= 0` or `Mdot_halo <= 0` (unformed/non-accreting progenitor) forces `S_uv = 0` rather than evaluating the formula (which would otherwise divide by zero).

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| UVB_enabled | `reionization.UVB_enabled` | True | Master on/off switch |
| z_reion | `reionization.z_reion` | 7 | Reionization redshift |
| gamma | `reionization.gamma` | 15 | Sharpness of the ionization transition |
| omega | `reionization.omega` | 2 | Power-law index of the suppression functional form |

## Merger tree input ([`utils.py`](ashvini/utils.py), [`pymctrees_adapter.py`](ashvini/pymctrees_adapter.py))

Every equation above consumes a halo's mass/growth-rate history along cosmic time, from either an HDF5 file (`basics.tree_source: file`) or a live [pymctrees](https://github.com/doctorcbpower/pymctrees) forest (`basics.tree_source: pymctrees`). See [README.md](README.md#input-halo-merger-trees) for the config surface. Notable conversions in the adapter: pymctrees' Msun/h masses and backward-in-time ordering are converted to plain Msun / forward-chronological; growth rates are finite-differenced against Ashvini's own (Planck18) `time_at_z`, not pymctrees' internal cosmology; and pre-formation placeholder masses (below pymctrees' resolution) are zeroed to represent "halo does not exist yet," which is what the `M_halo <= 0` guards throughout this document are protecting against.

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| tree_source | `basics.tree_source` | `"file"` | `"file"` reads `tree_file`; `"pymctrees"` generates a forest live each run |
| tree_file | `basics.tree_file` | `./data/inputs/merger_trees.h5` | HDF5 merger-tree input (used when `tree_source: file`) |
| mass_bin | `basics.mass_bin` | 1.0e10 Msun | z=5 halo mass bin label, used for output naming and `pymctrees` live generation's default `m_res` |
| dir_out | `basics.dir_out` | `./data/outputs/` | Output directory |

## Numerical scheme

`run1()`/`run_forest()` (the production path) integrate the same five ODEs as `run1_scalar()` (the `solve_ivp` reference), but using closed-form/quadrature updates chosen per-equation:

- **Stellar mass and stellar metals**: not self-referential within a step (pure integrals of quantities frozen at the step start) -> exact 5-point Gauss-Legendre quadrature of `1/t_ff(z(t))` over the step.
- **Gas mass, gas metals, dust mass**: each is affine, `dy/dt = forcing - decay*y`, with forcing/decay frozen at the step midpoint redshift -> solved exactly via an exponential integrator (`_linear_ode_step`).
- **Black hole mass**: piecewise (gas-supply-limited vs. Eddington-limited, with a possible regime crossover mid-step) -> solved exactly per-regime rather than approximated.

This lets the integration vectorise across every halo in a forest simultaneously (they share a redshift grid), rather than looping per-halo with `solve_ivp`. `tests/test_run1.py` cross-validates the vectorised path against the scalar reference.

## Citation

See [README.md](README.md#citation) for the papers to cite.
