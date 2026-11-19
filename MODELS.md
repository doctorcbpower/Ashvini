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

Gas-phase metal mass. All metallicities are absolute mass fractions inside the code; the config states them in units of Z_sun (`constants.Z_SUN = 0.015`) and the run header prints both:

```
dM_Z,gas/dt = Z_IGM*Mdot_acc - (eps_ff/t_ff) * M_Z,gas + Z_yield*SFR_wind
              - eta_w(z,M_halo,Z_star) * (M_Z,gas/M_gas) * SFR_wind        (last term omitted if M_gas = 0)
```

Stellar metal mass (new stars inherit the gas-phase metal fraction at formation, `SFR * M_Z,gas / M_gas = (eps_ff/t_ff) * M_Z,gas`):

```
dM_Z,star/dt = (eps_ff/t_ff) * M_Z,gas
```

`SFR_wind` is the same quantity as in the gas equation (zero, present or delayed, set by `supernova.type`). In the vectorised integrator the stellar mass and stellar metals formed in a step follow from mass balance over the step (for `dy/dt = F - D y`, the time integral of `D y` is `y0 + F dt - y1`), so gas + stars + removed mass balance to round-off and a reservoir that is empty at the start of a step forms stars in it.

Stellar metallicity `Z_star = M_Z,star/M_star` (0 if `M_star=0`) feeds back into the mass-loading factor and dust wind-loss below.

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| Z_IGM | `metallicity.Z_IGM_Zsun` | 1.0e-3 | Metallicity of inflowing IGM gas, in units of Z_sun (= 1.5e-5 as a mass fraction) |
| Z_yield | `metallicity.Z_yield` | 0.06 | Mass of metals returned per unit mass of stars formed (a mass ratio, not in Z_sun) |

## Dust ([`dust.py`](ashvini/dust.py))

```
growth_rate       = Y_d * SFR_wind
SNe_rate          = nu_SN * SFR_wind                         (supernovae per unit time)
dust_loading      = 1 - exp(-(M_gas/M_crit)^Alpha)           (saturates to 1 for M_gas >> M_crit)
destruction_rate  = M_swept * (M_dust/M_gas) * dust_loading * SNe_rate
wind_loss         = (M_dust/M_gas) * eta_w(z,M_halo,Z_star) * SFR_wind    (both loss terms 0 if M_gas <= 0)

dM_dust/dt = growth_rate - destruction_rate - wind_loss
```

`SFR_wind` is the same star formation rate that drives the supernova wind in the gas and metal equations: zero for `supernova.type: "no"`, the present rate for `"instantaneous"`, the rate `supernova.delay_time` ago for `"delayed"`. `nu_SN` is the number of core-collapse supernovae per Msun of stars formed, computed from `supernova.imf` over `supernova.m_sn_min`-`m_sn_max` (Salpeter 8-40 Msun: 6.8e-3; 8-100 Msun: 7.4e-3). The swept-up mass `M_swept` carries the dust-to-gas ratio, so each supernova destroys `M_swept * M_dust/M_gas`, not a fixed mass. (Before 2026-11-02 the rate was `Gamma * SFR_past / M_star,past`, which equals `Gamma / dt_step` and depended on the tree grid, and there was no `M_dust/M_gas` factor.)

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| M_swept | `dust.m_swept` | 1.0e3 Msun | ISM mass swept/destroyed per SN (`= epsilon_d * M_ISM`) |
| Y_d | `dust.dust_yield` | 0.004 | Dust mass produced per unit (delayed) SFR |
| nu_SN | `supernova.imf`, `m_sn_min`, `m_sn_max` | salpeter, 8, 40 Msun | Sets the supernova rate per Msun formed (the upper limit is a parameter; it does not change the mass-loading factor) |
| Alpha | `dust.dust_alpha` | 8 | Steepness of the destruction-efficiency saturation curve |
| M_crit | `dust.m_crit` | 1.0e5 Msun | Gas-mass scale at which destruction efficiency saturates |

## Supernova feedback ([`supernovae_feedback.py`](ashvini/supernovae_feedback.py))

Metallicity modulation (smooth logistic step from strong feedback at low `Z_star` to weak feedback at high `Z_star`; the step center/width/floor/ceiling below are hard-coded, not exposed in the YAML; Z_star is converted to units of Z_sun first, so the step is centred on 0.1 Z_sun):

```
f_Z(Z_star) = sigmoid((Z_star/Z_sun - 0.1)/0.01) * (0.25 - 1) + 1
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
| delay_distribution | `supernova.delay_distribution` | `"single_delay"` | The distribution of times between a star forming and exploding: `"single_delay"` (every star after `delay_time`) or `"stellar_lifetimes"` (each massive star after its own lifetime) |
| wind_model | `supernova.wind_model` | `"ashvini"` | The wind model, i.e. how much gas each supernova removes: `"ashvini"` (the mass loading above), `"momentum"` (momentum-driven) or `"energy"` (energy-driven) |

### Delay between star formation and explosion, and the wind model (`components/delay.py`, `components/wind_model.py`)

**Delay distribution `"stellar_lifetimes"`.** Each massive star that explodes as a supernova, with initial mass `m` drawn from the IMF between `m_sn_min` and `m_sn_max` (`n_progenitor_masses` bins of equal probability), explodes `t_MS(m, Z)` after it forms, where `t_MS` is its main-sequence lifetime (`stellar_lifetimes.py`; Georgy et al 2013 at `lifetime_Z` = 0.002, Ekstrom et al 2012 at 0.014, interpolated in log Z between them; rotating models; the lifetime ends when the central hydrogen fraction falls to 1e-4). Let `F(s)` be the fraction of these stars that have exploded a time `s` after forming, and `G` the integral of `F`. Then the fraction of the stars formed (at a constant rate) in step `k` whose supernovae go off in step `j` is `K[j,k] = [G(T_j - T_{k-1}) - G(T_j - T_k) - G(T_{j-1} - T_{k-1}) + G(T_{j-1} - T_k)] / dT_k`, where `T` are the cosmic times at the step boundaries and `dT_k` is the length of step `k`; this is exact when the star formation rate is constant within a step. The star formation rate that drives the wind in step `j` is `sum_k K[j,k] S_k / dt_j` over earlier steps, with `S_k` the mass of stars formed in step `k`. Supernovae that would go off within the step in which their stars formed are counted at the start of the next step. The same rate drives the loss of gas, metals and dust. With a single lifetime this is the single delay `delay_time` (to within one step). Used only when `supernova.type` is `"delayed"`.

**Wind models.** The mass loading `eta_w` is the gas removed per unit mass of stars formed, `eta_w = nu_SN * (gas removed by one supernova)`, with `nu_SN` the number of supernova progenitors per Msun of stars formed for the chosen IMF and mass range:

```
ashvini   eta_w from the calibrated formula above
momentum  eta_w = nu_SN * f_p * p_SN / v_esc                 (p_SN = 3e5 Msun km/s, the momentum one supernova delivers)
energy    eta_w = nu_SN * 2 * eps_fb * E_SN / v_esc^2        (E_SN = 1e51 erg, the energy of one supernova;
                                                              eps_fb = the fraction of it that unbinds gas)
```

`v_esc = sqrt(2) V_c`, with `V_c` the circular speed at the virial radius of a halo whose mean density is 200 times the critical density (`halo.py`). The Ashvini mass loading falls as `v_esc^-1`, like a momentum-driven wind: it equals the momentum-driven wind with `f_p` = 0.7646 at zero stellar metallicity to 0.3 percent over 1e6-1e10 Msun and z = 5-15 (for 8-100 Msun progenitors). `presets.stochocc_smooth_control()` sets this up together with Stochastic_Occupation's star formation rate. The grid integrator does not limit the gas removed to the gas present (a negative gas mass is set to zero afterwards); the event integrator does.

### Event integrator (`integrators/events.py`, `main.run_events`)

Follows individual supernovae for one halo at a time. Between supernovae the gas, stars and metals are advanced exactly (`dM/dt = S - kM`, solved in closed form, with `k` the star formation rate per unit gas mass at the start of the substep and the gas supply `S = f_b max(dM_h, 0)` reduced by the equilibrium UV suppression of Okamoto et al 2008). The stars formed in a substep are passed to `samplers.py`, which returns the supernovae they will produce. Each supernova acts at the moment it explodes, using the halo mass and escape speed at that moment: it adds its metals and dust (`Z_yield` and `dust_yield` per Msun of stars, divided by `nu_SN` per supernova), destroys dust (`weight * m_swept * M_d/M_g * f_gas`, at most the dust present), and removes `min(M_gas, weight * (gas removed by one supernova))` of gas, with metals and dust in proportion to the gas. Here `weight` is the fraction of a whole supernova that the event stands for. Two ways of choosing the supernovae: `"smooth"` (the expected supernovae, as fractions of a supernova, on K = 100 representative masses, released in batches of 1 Myr; no randomness) and `"stochastic"` (a Poisson number of whole supernovae with random masses). Not included: black holes. The gas supply and the UV suppression use the equilibrium form, unlike the grid integrator's form that depends on the growth rate.

#### Optional physics of the event integrator (Config section `events`, all off by default)

Each is compared with Stochastic_Occupation's smooth model (`tests/test_event_options.py`; reference numbers in `tests/reference_outputs/stochocc_options.json`). `tau` below is the Config value (a number of crossing times) times `r_vir / v_esc` at the time named, where `r_vir / v_esc` is the time to travel the virial radius `r_vir` at the escape speed `v_esc`.

| Key | What it does |
|---|---|
| `cooling_T_vir_K`, `cooling_mass`, `reionised_below_z`, `cooling_T_vir_reionised_K` | Cooling threshold: gas forms stars only when the halo's virial temperature `T_vir` (mu = 1.22) is above `cooling_T_vir_K`, or its mass above `cooling_mass` (which replaces the temperature condition); below redshift `reionised_below_z` it must also exceed `cooling_T_vir_reionised_K` (gas in smaller haloes cannot cool after reionisation). Tested at the start of each substep. |
| `inflow_delay_in_crossing_times` | Gas takes a time `tau` (at the start of the substep) to reach the galaxy: the supply passes through a reservoir that empties at rate `1/tau`, solved exactly. |
| `recycling_fraction`, `recycling_delay_in_crossing_times`, `recycling_fraction_Z` | This fraction of the gas a supernova removes falls back into the galaxy after `tau` (at the supernova); `recycling_fraction_Z` (default: the same) is the fraction of the removed metals and dust that fall back. |
| `stored_energy_loss_time_myr` | Energy kept for later: when a supernova finds less gas than its energy could unbind, the unused energy `(wanted - removed)/2 * v_esc^2` is stored, decays over `stored_energy_loss_time_myr` (`inf` = never), and removes gas as it arrives, before it reaches the galaxy. Needs the energy-driven wind. |
| `reachable_gas_fraction` | One supernova can reach only this fraction of the gas present. Because the limit applies to each supernova separately, the answer depends on how the expected supernovae are divided into separate explosions (the smooth sampler's `K`). |
| `coupling_change_mass`, `eps_fb_above_change_mass` | Above this halo mass at the time of the explosion the energy-driven wind uses a different coupling `eps_fb_above_change_mass` (supernovae confined in massive minihaloes). Needs the energy-driven wind. |
| `wind_delay_in_crossing_times` | The gas a supernova will remove is fixed when it explodes, but leaves the galaxy a time `tau` later (limited to the gas then present); metals and dust are lost with it then. The gas can form stars in the meantime. |

Stochastic_Occupation does not allow several of these together with chemistry (or with the momentum-driven wind). The event integrator always follows metals and dust and allows every combination, so those combinations are not compared with it.

## Disc size and a star formation threshold on surface density (`components/disc.py`, Config section `disc`, off by default)

With the option off the gas is one bulk reservoir. With `disc.enabled: True` the gas is taken to form an exponential disc, `Sigma(r) = Sigma_0 exp(-r / R_d)`, holding all of the gas mass `M_gas`, so `Sigma_0 = M_gas / (2 pi R_d^2)`. The scale length follows Mo, Mao & White (1998) without their corrections for the halo concentration and for the disc's pull on the halo:

```
R_d = (j_d / m_d) * lambda * R_vir / sqrt(2)
```

where `lambda` is the halo spin parameter (`spin_parameter`, 0.035, constant in time), `j_d / m_d` is the disc's share of the halo's angular momentum divided by its share of the mass (`angular_momentum_ratio`, 1) and `R_vir` is the virial radius (mean density 200 times critical). For `j_d / m_d = 1` this is the disc size used by the reservoir model, and the two virial radii agree to 0.3 percent.

Gas forms stars only where the surface density is above `critical_surface_density_Msun_pc2` (`Sigma_crit`, default 5 Msun/pc^2; 0 = no threshold). For an exponential disc that is every radius inside `R_crit = R_d ln(Sigma_0 / Sigma_crit)`, which holds the fraction

```
f = 1 - (1 + x) exp(-x),    x = ln(Sigma_0 / Sigma_crit) = R_crit / R_d      (f = 0 if Sigma_0 <= Sigma_crit)
```

of the gas. The star formation rate is the usual rate applied to that gas only: `SFR = f * eps_ff / t_ff * M_gas` (grid integrator: `f` from the gas and halo mass at the start of the step; event integrator: at the start of the substep). `f` rises smoothly from 0 as the disc becomes dense, so star formation does not switch on and off from step to step. `components.disc.diagnostics` gives `R_d`, `Sigma_0` and `f` for a finished run.

**Star formation law from the disc (`disc.star_formation_law`).** `"halo_free_fall"` (the default) is the usual rate on the gas above the threshold, as above. `"schmidt_kennicutt"` sets the rate from the disc: `Sigma_SFR = A Sigma_gas^N` (Kennicutt 1998; `schmidt_kennicutt_normalisation` A = 2.5e-4 Msun/yr/kpc^2 for `Sigma_gas` in Msun/pc^2, `schmidt_kennicutt_slope` N = 1.4), integrated over the exponential disc, `SFR = 2 pi A Sigma_0^N R_d^2 / N^2`. It has no threshold. The relation is calibrated on present-day galaxies and here applied at every redshift and mass. `components.disc.implied_efficiency` gives the efficiency per halo free-fall time that a disc law amounts to at the present halo and gas mass.

A rate `eps_d * M_gas / t_disc`, with `t_disc = R_d / V_c` the disc's dynamical time, is deliberately not an option. For `R_d` proportional to `R_vir` and `V_c` the halo circular speed, `R_vir / V_c = 1 / (10 H(z))`, so `t_disc` is a fixed multiple of `1 / H(z)` and the law is exactly the usual rate with efficiency `eps = eps_d * 0.141 * 10 sqrt(2) / ((j_d / m_d) * spin)`: 0.855 for `eps_d` = 0.015 and spin 0.035, for every halo and redshift. It would differ only if the spin varied between haloes or the circular speed at `R_d` differed from the halo's. `tests/test_disc.py` checks this.

**Can the usual law be made consistent by changing its efficiency?** For the disc dynamical time law, yes: it is the same law (above). For the Schmidt-Kennicutt law, not with one number (`scripts/disc_star_formation_laws.py`, three haloes per mass; output in `scripts/output/disc_star_formation_laws.txt`). The efficiency per halo free-fall time that it amounts to is about 0.06 at 1e8 Msun, 0.12 at 1e10 and 0.16 at 1e11 (median over time), so it rises with mass, and within each history it changes by a factor of 4 to 7 (5th to 95th percentile; at 1e10 from 0.02 at z = 20 to 0.4 at z = 6) because it grows with the gas surface density (`eps` proportional to `Sigma_gas^(N-1)`). A single efficiency chosen to give the same final stellar mass is 0.12 to 0.14 at 1e8, 0.1 to 1.1 at 1e10 (the stellar mass is not a monotonic function of the efficiency there, so the match is not unique) and 0.5 to 0.7 at 1e11, i.e. different for each mass. With it the gas mass at z = 5 is 0.6 to 0.7 of the Schmidt-Kennicutt value at 1e8, 0.3 to 3 at 1e10 and 0.7 to 1.2 at 1e11, and the star formation rate at early times differs by factors of 3 to several hundred. The gas metallicity is matched to 1 to 30 percent. A constant efficiency therefore reproduces the final stellar mass at a given mass but not the gas content or the history. An efficiency that depends on the gas surface density (and so on mass and time) is needed.

**Time step.** A star formation rate much faster than the supernova delay or the step makes the grid result depend on the redshift grid (feedback acts a supernova delay late relative to the step). `GridIntegrator` warns when star formation consumes more than about 26 percent of the gas in a step. The Schmidt-Kennicutt law is first order in the step: 3 percent below the converged value at the resolution of the test trees. The event integrator shortens its substeps to a tenth of the star formation time when the disc option is on.

What this leaves out: the concentration and disc-gravity corrections (the fuller formula in the DiscModels package gives a disc 0.55-0.70 times as large for the cases tried, which `angular_momentum_ratio` can mimic), a spin that changes with time or differs between haloes, the disc's stability (Toomre) condition, molecular hydrogen, and the disc's effect on the dark matter halo. The disc is a test of what a disc-based star formation does, not a disc galaxy model. `scripts/disc_threshold_subset.py` (the threshold; it changes the stellar mass by 3 percent or less for thresholds up to 10 Msun/pc^2) and `scripts/disc_star_formation_laws.py` (the rate law; the Schmidt-Kennicutt law gives 1.5 to 11 times the stars of the usual law at the fiducial parameters, 1e8 to 1e11 Msun) compare them with the usual law on a few haloes.

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
| Z_gas_max | `black_holes.seeding.pop3.Z_gas_max_Zsun` | 1.0e-4 | Near-primordial gas required (units of Z_sun) |
| M_seed | `black_holes.seeding.pop3.M_seed` | 1.0e2 Msun | Light seed (Pop III remnant) |
| enabled | `black_holes.seeding.direct_collapse.enabled` | True | Direct-collapse channel on/off |
| z_min | `black_holes.seeding.direct_collapse.z_min` | 10 | Redshift floor |
| M_halo_min | `black_holes.seeding.direct_collapse.M_halo_min` | 1.0e7 Msun | Atomic-cooling halo threshold |
| Z_gas_max | `black_holes.seeding.direct_collapse.Z_gas_max_Zsun` | 1.0e-5 | Pristine gas required to suppress H2 fragmentation (units of Z_sun) |
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
| growth_cap_enabled | `black_holes.sigma_feedback.growth_cap_enabled` | False | Hard `M_BH` growth cap, see below. Requires `enabled: True` too. |

**Important caveat about `enabled` alone, checked directly (see `tests/test_sigma_feedback.py::test_sigma_feedback_wind_alone_does_not_cap_bh_growth_near_m_sigma`): by itself, this does not hold `M_BH` near `M_sigma`.** "Self-regulation" in King (2003, 2005)/PZNK11 refers to the BH's own growth eventually being throttled once the wind is strong enough to expel/starve the gas supply feeding it -- but `enabled: True` alone only adds a gas-removal *wind* term (`Mdot_wind = eta_agn_eff * Mdot_BH`), not a direct cap on `M_BH` itself. Regulation, if any, can only happen indirectly through the wind depleting `gas_mass` enough to lower the gas-supply-limited growth rate (`epsilon_BH/t_ff * gas_mass`) that feeds `M_BH`. Checked over a long (z=25 -> z=0.5), sustained-accretion integration: that indirect loop is nowhere near strong enough -- cosmological gas resupply vastly outpaces even a maximal wind (`eta_agn=1.0`), so `M_BH` keeps growing over two orders of magnitude past `M_sigma`, not settling near it. `growth_cap_enabled` below (off by default) is what actually closes this gap.

#### Hard growth cap (`growth_cap_enabled`)

PZNK11 do not give an equation for *how* growth halts once `M_BH` approaches `M_sigma` -- in their own words, **"precisely how this last step occurs is still the most unclear part of any theory."** What they do give (eq. 21-22) is a bound on the eventual overshoot, reached once the AGN outflow "presumably halts further accretion" after propagating through the bulge on its dynamical time:

```
ceiling(M_halo, z) = M_sigma * (1 + 0.41 * sigma_200 / h(z))
```

`sigma_200 = sigma(M_halo, z) / (200 km/s)`, `h(z)` the dimensionless Hubble parameter (`utils.h_of_z`). `growth_cap_enabled: True` (`black_holes_growth.growth_ceiling`) applies this as a **hard cutoff**: `dM_BH/dt = 0` once `M_BH >= ceiling`, implemented by clipping each step's growth-update result to the ceiling (exact, not an approximation -- `M_BH` is monotonically non-decreasing in this model, so clipping the endpoint is equivalent to freezing the ODE the instant it would cross the ceiling). Because `ceiling` grows with the halo (via `sigma`), a previously-capped BH can resume growing later if its halo grows enough to raise the ceiling -- this is a deliberate consequence of tying the cap to the halo's own evolving properties, not a bug.

**The hard cutoff itself is a disclosed modelling choice, not an equation read directly off the paper** -- PZNK11 gives the ceiling value, not the transition's functional form. Checked directly (`tests/test_sigma_feedback.py::test_growth_cap_enabled_actually_holds_bh_mass_near_ceiling`): under the same sustained-accretion, maximal-wind conditions where the wind-only mechanism let `M_BH` grow two orders of magnitude past `M_sigma`, the hard cap holds `M_BH` at >90% of its ceiling throughout.

Out of scope for this mechanism (potential follow-ups): an NFW-profile version of `sigma`/`M_sigma`; a smooth (rather than hard-cutoff) transition at the ceiling; and the competing nuclear-cluster feedback channel (Nayakshin, Wilkinson & King 2009).

## Alternative BH growth model: `growth_model: "hobbs_slimdisk"` ([`black_holes_growth_slimdisk.py`](ashvini/black_holes_growth_slimdisk.py), [`spin.py`](ashvini/spin.py))

Selected via `black_holes.growth_model` (default `"pznk11_freefall"`, the
model above -- unchanged). Ported from the 2026 "Differential Growth"
paper's standalone prototype; see
[`docs/2026_paper_session_code_catalogue.md`](docs/2026_paper_session_code_catalogue.md)
for provenance and what was deliberately *not* ported (a stochastic,
time-varying compaction boost -- see below). Genuinely different physics
from the default model above, not a bug-fix relationship -- both are
physically motivated but make different choices; see the catalogue doc's
comparison table.

**Nuclear accretion supply** (Hobbs, Power, Nayakshin & King 2012): a
free-fall estimate on the *total enclosed mass*, not naive Bondi-Hoyle
(`Mdot ~ M_BH^2`), which is the wrong estimator whenever `M_BH` does not
dominate the local enclosed mass -- exactly the light-seed case:

```
t_ff(M_enc) = sqrt(R_nuc^3 / (2*G*M_enc)),   M_enc = M_BH + M_gas + M_star
Mdot_acc = eta_acc * (M_gas / t_ff(M_enc)) * Phi_hat^2
```

`Phi_hat` (`black_holes.slimdisk.compaction_boost`) is a Booth & Schaye
(2009)-style density boost, anchored to Lapiner, Dekel & Dubois (2021)
compaction events. **Simplification versus the paper's own analysis
scripts**: the paper explored a stochastic, time-varying `Phi_hat(t)`
(an Ornstein-Uhlenbeck process); that is not reproduced here, since it
doesn't fit this package's "coefficients frozen over a step" closed-form
integration scheme (see [Numerical scheme](#numerical-scheme)) --
`compaction_boost` here is a single fixed multiplier for a given run.

**Radiative efficiency / spin**: exposed as a genuine free parameter
(`black_holes.slimdisk.a_star`, BH spin, mapped to `epsilon` via the
Novikov-Thorne ISCO relation, `spin.epsilon_from_spin`) -- unlike the
default model above, whose Eddington-rate constant has no `epsilon`
dependence at all.

```
kappa_edd = [(1-epsilon)/epsilon] * 4*pi*G*m_p / (sigma_thomson*c)   [1/Gyr]
Mdot_Edd,std(M_BH) = kappa_edd * M_BH
```

**Graded super-Eddington cap** (Watarai et al. 2000; Madau, Haardt &
Dotti 2014; Lupi et al. 2024 slim-disc picture) -- a threshold ratio
`r_crit`, not a hard multiplier applied at every mass:

```
cap(M_BH) = Mdot_Edd,std(M_BH)   if Mdot_acc/Mdot_Edd,std(M_BH) <= r_crit
            Mdot_acc / r_crit    otherwise

dM_BH/dt = min(Mdot_acc, cap(M_BH))
```

With `Mdot_acc` frozen over a step (same pattern as every other
time-varying coefficient in this package), this is a three-regime
piecewise ODE in `M_BH` -- constant (slim-disc-throttled, small `M_BH`),
exponential (standard-Eddington-limited, intermediate `M_BH`), constant
(gas-supply-limited, large `M_BH`) -- solved exactly in closed form by
`bh_growth_step_slimdisk` (crossing times between regimes found
analytically, not approximated), verified against a numerical
`solve_ivp` reference and against the default model's own
`_bh_growth_step` in the limit `r_crit -> infinity` (see
`tests/test_black_holes_growth_slimdisk.py`).

**AGN feedback**: King (2003) energy-driven wind, a direct, constant-`epsilon_f`
coupling (not the default model's `M_BH/M_sigma` logistic switch):

```
Mdot_wind = (2*epsilon_f*c^2 / sigma^2) * Mdot_BH
```

`sigma` reuses the default model's isothermal-sphere `velocity_dispersion(M_halo, z)`.

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| growth_model | `black_holes.growth_model` | `"pznk11_freefall"` | `"hobbs_slimdisk"` selects this model instead |
| eta_acc | `black_holes.slimdisk.eta_acc` | 0.01 | Nuclear accretion efficiency |
| R_nuc_pc | `black_holes.slimdisk.R_nuc_pc` | 100.0 | Fixed nuclear radius for the free-fall estimate, pc |
| compaction_boost | `black_holes.slimdisk.compaction_boost` | 1.0 | Fixed `Phi_hat` multiplier (>=1); no time-varying version, see above |
| r_crit | `black_holes.slimdisk.r_crit` | 8.0 | Supply/standard-Eddington ratio threshold for the slim-disc cap |
| epsilon_f | `black_holes.slimdisk.epsilon_f` | 5.0e-4 | King (2003) AGN wind coupling |
| a_star | `black_holes.slimdisk.a_star` | 0.5 | BH spin -> radiative efficiency via `spin.epsilon_from_spin` |

A companion utility, [`seed_mass_function.py`](ashvini/seed_mass_function.py),
converts a critical-seed-mass boundary curve (`M_seed,crit(M_halo)`, the
seed mass below which strict-Eddington growth under this model cannot
reach a target `f_BH` by some anchor redshift -- see the paper) into
predicted boosted-fractions for standard seed-formation channels
(Pop III / runaway-collision / direct-collapse), given illustrative
log-normal-in-log10(mass) population models for each. Not a rigorously
derived population synthesis -- see the module docstring.

## Critical seed mass ([`critical_seed.py`](ashvini/critical_seed.py))

General-purpose companion to `seed_mass_function.py` above: given a halo's
assembly history, finds the black-hole seed mass `M_seed,crit(M_halo)` at
which strict Eddington-limited growth (`growth_model: "hobbs_slimdisk"`,
super-Eddington cap disabled -- see below) reaches a target
overmassiveness by a chosen anchor redshift,

```
M_BH(z_anchor; M_seed,crit) = f_BH * M_star(z_anchor),   f_BH = 0.5 fiducial
```

`M_seed < M_seed,crit` means Eddington-limited growth alone cannot reach
the target by `z_anchor` (a super-Eddington episode, earlier formation, or
a lower effective radiative efficiency would be needed to explain an
observed overmassive black hole); `M_seed >= M_seed,crit` means it can.
This is the light-seed-problem boundary itself, computed directly against
`main.run_forest()` rather than the paper-specific reservoir model of
`paper_reservoir.critical_seed_paper` -- any caller with a halo assembly
history (`halo_mass`, `halo_mass_rate`, `redshift`, the same triple
`build_forest_live`/`build_forest_for_bin` return) can use it, not only
the 2026 paper's own bespoke gas-supply model.

**Strict Eddington-limited growth means `r_crit -> infinity`, not
`r_crit=1`.** `black_holes_growth_slimdisk.bh_growth_step_slimdisk`'s
closed-form update has three regimes with boundaries
`M1=A_bh/(kappa_edd*r_crit)` and `M2=A_bh/kappa_edd=M1*r_crit`; at
`r_crit=1`, `M1=M2` exactly, so the genuine Eddington-limited
(exponential, `kappa_edd`-dependent) regime between them has zero width
and growth equals the raw, uncapped nuclear supply rate `A_bh` at every
black hole mass -- i.e. `r_crit=1` *removes* the Eddington cap entirely,
the opposite of what an earlier version of this module assumed. This was
a real bug (found and fixed 2026-09-15, in both this module and the
sibling `paper_reservoir.critical_seed_paper`, which reimplements the
identical bisection against `paper_reservoir.run_reservoir_paper` instead
of `main.run_forest`): `critical_seed()`'s `r_crit` parameter used to
default to `1.0`, silently computing an *uncapped*-growth seed-mass
boundary while claiming to compute the strict-Eddington one.
`STRICT_EDDINGTON_R_CRIT = 1e12` is now the default, verified directly
against `bh_growth_step_slimdisk` to reproduce pure Eddington-exponential
growth (no crossover into the supply-limited regime) over the seed-mass
range this search brackets. Passing a smaller, finite `r_crit` computes a
*different*, still well-defined quantity -- a seed-mass boundary under a
permitted super-Eddington episode -- just not "M_seed,crit" in the
paper's strict-Eddington sense.

Implementation: vectorised bisection in log-seed-mass, one root per halo,
against `excess(seed_mass) = M_BH(z_anchor) - f_BH*M_star(z_anchor)`.
Well-posed because growth is monotonically non-decreasing in seed mass at
`r_crit=STRICT_EDDINGTON_R_CRIT` (the two-regime closed form is monotonic
in its initial condition), provided AGN feedback's own dependence on seed
mass (via `epsilon_f`) doesn't overwhelm that monotonicity -- true for any
realistic `epsilon_f`, since the AGN wind is a small perturbation on the
gas budget, not a dominant term (see AGN feedback above). Seeding is
injected at each halo's own first-resolved step (not a fixed halo-mass
threshold) by temporarily repointing the existing
`black_holes.seeding.halo_mass_threshold` channel at a per-halo trial
seed-mass array and restoring the original config afterward -- no core
seeding-channel changes were needed for this.

| Parameter | Meaning | Default |
|---|---|---|
| f_bh | Target overmassiveness `M_BH(z_anchor)/M_star(z_anchor)` | 0.5 |
| seed_mass_lo, seed_mass_hi | Bisection bracket, Msun | 1.0, 1.0e8 |
| n_iter | Bisection iterations | 50 |
| m_halo_min | Halo-mass floor for the seeding-injection trick above, Msun | 1.0 |
| r_crit | Super-Eddington cap threshold; `STRICT_EDDINGTON_R_CRIT` (1e12) for the paper's strict-Eddington sense | `STRICT_EDDINGTON_R_CRIT` |

Returns `M_seed_crit` (NaN where the bracket doesn't resolve a root: either
`seed_mass_hi` never reaches the target, or `seed_mass_lo` already exceeds
it -- both returned as separate boolean masks so a caller can tell which
failure mode fired and widen the bracket accordingly) plus `**overrides`
(`a_star`, `epsilon_f`, `eta_acc`) forwarded unchanged to `run_forest()`.
Not yet re-exercised against real data since the `r_crit` fix, unlike the
sibling `paper_reservoir.critical_seed_paper`, which was stress-tested
extensively as part of the 2026 paper session (see
`docs/2026_paper_session_code_catalogue.md`) -- callers relying on this
module should re-verify their own results were computed after the fix.

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

**Audit note (2026-09-20, [docs/UV_SUPPRESSION_AUDIT.md](docs/UV_SUPPRESSION_AUDIT.md)).** The correction term equals the rate at which the equilibrium gas fraction `s(M/M_c) f_b` falls as `M_c(z)` grows, and inflow cannot go negative, so the term stops accretion rather than removing gas. Measured on Zhang-Hui trees it lets essentially no accretion through after z = 5 for `M_halo < 1e10` and about half at 1e11, and it makes the z = 0 stellar mass depend on the tree time step. The `epsilon(z)` in the code equals `-d ln M_c / d z` for z <= 5 but not above (0.09 of the true value at z = 7, zero for z >= 8); the formula written for `epsilon` above does not match either. `paper_reservoir.uv_suppression_mass_only` (the `s(mu, omega)` factor alone) is the form the MVM uses.

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| UVB_enabled | `reionization.UVB_enabled` | True | Master on/off switch |
| z_reion | `reionization.z_reion` | 7 | Reionization redshift |
| gamma | `reionization.gamma` | 15 | Sharpness of the ionization transition |
| omega | `reionization.omega` | 2 | Power-law index of the suppression functional form |

`reionization.form` chooses how the suppression is applied to the gas accretion. `"rate_dependent"` (default) is the form above, with the term that depends on the halo growth rate; `docs/UV_SUPPRESSION_AUDIT.md` found that it stops essentially all accretion after z = 5 in haloes below about 1e10 Msun and that it depends on the tree time step. `"equilibrium"` uses only `s(M / M_c(z), omega)`, the suppression of the equilibrium baryon fraction at the halo mass and redshift (Okamoto, Gnedin & Frenk 2008; the form the reservoir model uses, and the form the event integrator always uses): accretion is reduced, not stopped, and the result does not depend on the growth rate or the time step.

## Merger tree input ([`utils.py`](ashvini/utils.py), [`pymctrees_adapter.py`](ashvini/pymctrees_adapter.py))

Every equation above consumes a halo's mass/growth-rate history along cosmic time, from either an HDF5 file (`basics.tree_source: file`) or a live [pymctrees](https://github.com/doctorcbpower/pymctrees) forest (`basics.tree_source: pymctrees`). See [README.md](README.md#input-halo-merger-trees) for the config surface. Notable conversions in the adapter: pymctrees' Msun/h masses and backward-in-time ordering are converted to plain Msun / forward-chronological; growth rates are finite-differenced against Ashvini's own (Planck18) `time_at_z`, not pymctrees' internal cosmology; and pre-formation placeholder masses (below pymctrees' resolution) are zeroed to represent "halo does not exist yet," which is what the `M_halo <= 0` guards throughout this document are protecting against.

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| tree_source | `basics.tree_source` | `"file"` | `"file"` reads `tree_file`; `"pymctrees"` generates a forest live each run |
| tree_file | `basics.tree_file` | `./data/inputs/merger_trees.h5` | HDF5 merger-tree input (used when `tree_source: file`) |
| mass_bin | `basics.mass_bin` | 1.0e10 Msun | z=5 halo mass bin label, used for output naming and `pymctrees` live generation's default `m_res` |
| dir_out | `basics.dir_out` | `./data/outputs/` | Output directory |

## Minimal Viable Model of the galaxy and nuclear reservoirs ([`reservoir_stock.py`](ashvini/reservoir_stock.py))

The model behind the 2026 paper "Critical seed masses for massive black holes in the early Universe". It is separate from the ODE model above (`main.py`) and is integrated along a halo's main-progenitor history on a uniform cosmic-time grid with an explicit, causally ordered step. **Frozen**; the production calculation and its provenance are in [`docs/PRODUCTION_PROVENANCE.md`](docs/PRODUCTION_PROVENANCE.md). The module docstring is the authoritative statement of what is and is not included.

State: a black hole, two gas reservoirs (galaxy and nuclear), stellar mass at both scales, and a ledger of gas that fails the galaxy criterion. Defaults live in `paper_reservoir_params.yaml`.

* **Halo supply.** `Mdot_in = zeta_UV(M_halo, z) zeta_ch(M_halo) f_b Mdot_halo`, with the Okamoto et al. (2008) step for `zeta_UV` (`z <= 10`) and `zeta_ch = 1 - s(M_halo/M_hot, phi)`, `M_hot = 4e11 Msun`. `Mdot_halo` is a spline-derived rate from the tree history (`paper_reservoir.grumpy_halo_growth_rate`).
* **Delivery.** Each step's inflow is a cohort with `ln j ~ N(ln j_med, sigma_j)`, `j_med = lambda sqrt(2) V_vir R_vir`, `sigma_j = 0.5`. Galaxy cut (memoryless): the fraction with `j <= j_thresh(R_gal)`, `j_thresh(R) = sqrt(G M_enc R)`, joins the galaxy; the rest goes to the ledger. Nuclear transfer (with history): each cohort keeps the running maximum of `j_thresh(R_nuc)`, `R_nuc = 100 pc` fixed, `M_enc = M_BH + M_gas,nuc + M_star,nuc` (no dark matter); a later fall of the threshold moves nothing, a later rise releases gas. `transfer_nuc` is the mass moved per step.
* **Sinks.** `Mdot_star = eps_sf M_gas / t_ff` at both scales (`eps_sf = 0.015`, `t_ff = sqrt(R^3 / 2 G M_enc)`); `Mdot_BH = min(eta_acc M_gas,nuc / t_ff, Mdot_Edd)` with `eta_acc = 0.005` and an exact strict-Eddington step (`strict_eddington_step`). Both nuclear sinks draw on the same gas and are scaled together if their demand exceeds it.
* **Feedback.** Stellar outflow `eta_w(M_halo, z) Mdot_star,tot` (`eta_w = 5 (1e11.5/M_halo)^(1/3) (9/(1+z))^(1/2)`, metal-poor limit); momentum-driven AGN outflow `f_mom [eps/(1-eps)] (c/sigma_h) Mdot_BH`, `sigma_h` the isothermal halo velocity dispersion. Driven by the sinks actually achieved, removed from both gas reservoirs in proportion to gas mass, capped at the gas present, not recycled.
* **Critical seed.** `critical_seed_stock` finds `M_seed,crit` by bisection in `ln M_seed` on `F = M_BH(z=5)/(f_BH M_star,tot(z=5)) - 1`, so `M_seed,crit = f_BH M_star,tot / G_BH` with `G_BH = M_BH(z=5)/M_seed`.

Excluded: black-hole mergers, super-Eddington growth, nuclear dark matter, an evolving `R_nuc`, torques or angular-momentum loss, delayed feedback and wind recycling.

**Merger trees for the production run.** Zhang and Hui trees from `foraois` at `dz = 0.05`, `M_res = 1e4 Msun`. That timestep is outside the practical single-split compliance regime of the builder, so the trees are not timestep-converged (see the README and `foraois` `docs/PCH08_HIGH_Z_DIAGNOSTIC.md`); the effect on the black-hole results has not been quantified.

## Numerical scheme

`run1()`/`run_forest()` (the production path) integrate the same five ODEs as `run1_scalar()` (the `solve_ivp` reference), but using closed-form/quadrature updates chosen per-equation:

- **Stellar mass and stellar metals**: from mass balance over the step: for `dy/dt = F - D*y` the time integral of `D*y` is `y0 + F*dt - y1`, so the stars formed are `(B_sf / D) * (y0 + F*dt - y1)`, exact and conserving mass to round-off (`components/gas.py`, `components/chemistry.py`).
- **Gas mass, gas metals, dust mass**: each is affine, `dy/dt = forcing - decay*y`, with forcing/decay frozen at the step midpoint redshift -> solved exactly via an exponential integrator (`components/ode.linear_ode_step`).
- **Black hole mass**: piecewise (gas-supply-limited vs. Eddington-limited, with a possible regime crossover mid-step) -> solved exactly per-regime rather than approximated.

This lets the integration vectorise across every halo in a forest simultaneously (they share a redshift grid), rather than looping per-halo with `solve_ivp`. `tests/test_run1.py` cross-validates the vectorised path against the scalar reference.

## Citation

See [README.md](README.md#citation) for the papers to cite.
