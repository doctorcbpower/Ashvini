# Feeding-efficiency assumptions: what the paper can say (three points)

Working note, not manuscript text; the `.tex` is untouched. Scope, as clarified by the author: (1) why eps_sf = 0.015 is a reasonable nuclear star-formation efficiency, (2) what motivates eta_acc = 0.005 and how defensible a constant BH feeding efficiency is, (3) what the efficiency assumptions imply for the interpretation of M_BH/M_star, including that the 0.1 to 0.2 ratio is model-dependent. (`docs/mvm_discussion_alignment.md` was written against the manuscript's Discussion subsections, which was the wrong scope; its list of old-model conflicts is still valid but it is not this note.)

Sources of each statement are marked: **[model]** measured in the frozen MVM (claim-map ID); **[checked]** read by me in a source (see `docs/mvm_efficiency_literature_check.md`, fetch summaries "as extracted"); **[supplied]** stated in the author's message from the torque-model papers and not read by me; **[arithmetic]** computed by me from a supplied formula.

## 1. eps_sf = 0.015

What can be said.
* It is inherited from Menon, Balu & Power (2026), Table 1, where it multiplies a Hubble-scale clock tau_sf = 0.15 f_sf / H(z) **[checked]**. In the MVM the same number multiplies M_gas / t_ff(R) at each scale. The number was inherited and its meaning changed; at the nuclear scale (t_ff of order 1 Myr) the implied specific star-formation rate is about 100 times that of the galaxy-scale Hubble clock **[checked, arithmetic]**.
* Read as an efficiency per free-fall time, 0.015 lies in the observed range of about 0.01 to 0.1 (Krumholz & Tan 2007; McKee & Ostriker 2007; Utomo et al. 2018, from search results) **[checked]**. That literature does not constrain 100 to 300 pc nuclear regions.
* Its effect on the critical seed is modest: eps_sf x0.5 / x2 gives x0.77 / 1.20 (3e10) and x0.69 / 1.37 (3e11) **[model, R4]**, because the boundary tracks M_star,tot. Its effect on the accessible-regime ratio is direct: the no-feedback floor scales approximately as 1/eps_sf **[model, R8]**.

What cannot be said: that 0.015 is calibrated for nuclear star formation.

Possible external check **[supplied, unverified]**: Esquej et al. (2014), nuclear star formation on about 65 pc scales as a fraction of the 600 pc value, would test how much of the reservoir turns into stars. Not read; not used.

## 2. eta_acc = 0.005 and a constant feeding efficiency

What can be said.
* eta_acc is a new parameter of this paper. It was constrained only by the requirement eta_acc/eps_sf < f_BH (viability margin) and by the expectation that nuclear accretion is less efficient than star formation **[checked, session notes]**. Hobbs et al. (2012) and Hopkins & Quataert (2011) support the form of the supply (free-fall scaling; an angular-momentum bottleneck) but their abstracts give no value **[checked, abstracts only]**.
* Its effect on the critical seed is weak: eta_acc x0.1 / x10 gives x1.17 to 1.27 / x0.52 to 0.75 **[model, R6]**. Its effect on the accessible-regime ratio is strong (point 3).
* A constant eta_acc is an effective parameter standing for f_acc = [Mdot_torque(R0) / (M_g / t_ff)] x eps_m **[supplied]**: a resolved-variable inflow rate divided by the free-fall supply, times the fraction of inflow that survives to be accreted.

Indicative consistency check against a gravitational-torque inflow law **[supplied formula; arithmetic mine; not a calibration]**. Take Mdot_torque = eps_T f_d^{5/2} (M_BH/1e8)^{1/6} (M_d/1e9) (R0/100 pc)^{-3/2} (1 + f0/f_gas)^-1 Msun/yr, f0 = 0.31 f_d^2 (M_d/1e9)^{-1/3}, for a single reference nucleus with M_d = M_g = 1e8 Msun, f_d = f_gas = 1, M_BH = 1e6 Msun:

| eps_T | Mdot_torque (Msun/yr) at R0 = 200 pc / 100 pc | f_acc = Mdot_torque / (M_g / t_ff) | Mdot_torque / Mdot_Edd (eps = 0.1) at 200 pc |
|---|---|---|---|
| 0.05 | 0.0005 / 0.0014 | 2.1e-5 | 0.02 |
| 0.5 | 0.0049 / 0.0139 | 2.1e-4 | 0.22 |
| 5 | 0.049 / 0.139 | 2.1e-3 | 2.2 |

Recomputed: the equivalent f_acc is independent of R0 (both rates scale as R^-3/2), and reproduces the supplied values to within rounding. Mdot_Edd = 0.022 Msun/yr for 1e6 Msun.
* The MVM fiducial eta_acc = 0.005 is about a factor 2.4 above the eps_T = 5 value and about 24 above eps_T = 0.5, for this one nucleus. It is therefore at or beyond the generous end of the published normalisation range, not in the middle of it.
* f_acc from a torque law is not constant: it scales as f_d^{5/2} (M_BH)^{1/6} and falls with nuclear mass (recomputed: 4.5e-4, 2.1e-4, 8.4e-5 at M_g = 1e7, 1e8, 1e9, eps_T = 0.5). A constant eta_acc therefore overstates access in massive nuclei relative to this law. Direction only; not tested in the MVM.
* Limits: one system; f_d = 1 is the upper end (f_d^{5/2} lowers it); the calibration is for lower-redshift simulations, and fragmentation can compete in very gas-rich nuclei **[supplied]**; a nuclear supply fraction at z > 6 has no direct observational anchor that I know of.

What cannot be said: that eta_acc = 0.005 is calibrated, or that the torque model would reproduce or reverse any MVM result. Do not run an eps_T sweep in the MVM: eta_acc is not eps_T, and a sweep over the torque literature range would conflate two prescriptions (author's own instruction, 2026-09-20). The comparison above is a translation of one prescription into the other's units, not a sweep.

## 3. Consequences for the interpretation of M_BH/M_star

What can be said **[model]**.
* The ratio reached by light seeds in the accessible regime relaxes, in finite time, toward approximately eta_acc/eps_sf times a feedback factor (S3, S5, R8): no-feedback R/(eta_acc/eps_sf) = 0.6 to 0.99 (3e10), 0.3 to 0.94 (3e11), log-log slope 0.79 to 0.98 over a factor of about 80; with AGN and winds R/(eta_acc/eps_sf) falls to 0.08 to 0.44 (3e10, seed 1e3), slope about 0.56 to 0.57.
* The fiducial ratio of 0.1 to 0.2 is therefore a consequence of the adopted sink efficiencies and feedback, not an independent prediction of a characteristic BH-to-stellar mass scale (S3).
* The ratio is not a function of eta_acc/eps_sf alone: two grid points with equal ratio 0.667 give R = 0.546 and 0.401 (3e10) and 0.373 and 0.241 (3e11) (R8).
* Not at 3e13 (max/min across seeds 51) and not at the fiducial accessibility (no crossing; ratio spans 1e5 across seeds).

Regime correction for anything drawn from the torque framework. The supplied summary places the critical seed in the Eddington-limited regime and contrasts it with seed-independent, supply-limited torque growth. In the MVM this is not the case at the fiducial: the critical seed is supply-limited, with 0% Eddington-capped steps, 0.028 to 0.074 of about 20.8 e-folds used, and M_seed,crit within 3 to 8% of f_BH M_star,tot (H2, H3). Both the MVM and a torque law have a supply rate with weak dependence on M_BH, so the loss of seed memory in accessible systems does not require the constant-efficiency construction; what the construction fixes is the level of the ratio. The torque law is a useful contrast for how the supply depends on the nucleus (f_d, M_BH, M_d), not for a change of regime.

## 4. The three levels in this scope

1. Established by the model: given a nuclear supply and an effective feeding prescription, the critical seed is set by the stellar mass formed and by access to the nuclear reservoir (H2, H3, R5).
2. Established by the stress tests: the light-seed ratio in accessible systems follows from competing sinks and feedback, with a level set by the adopted efficiencies (S3, S5, R8); not a characteristic scale.
3. Unresolved: eta_acc has no independent calibration; eps_sf is inherited and applied per local free-fall time; a torque-regulated supply is one physically motivated replacement, and its use at z > 6 in gas-rich nuclei is uncertain. Whether replacing the constant with a torque law changes the critical-seed result qualitatively is a future question, not a result of this paper.

## 5. Status of the supplied literature material

Not read by me and not to be quoted until checked: eq. 2 and the eps_T range in Angles-Alcazar et al. (2017, arXiv:1603.08007); the eps_m factor of about 0.1 and the Hopkins et al. (2016) attribution; Simba (Dave et al. 2019); Diamond-Stanic & Rieke (2012); Esquej et al. (2014); Volonteri et al. (2015; the author list was not established); Hu et al. (2022, 2025) and the exponent p; Chen, Mo & Wang (arXiv:2509.03283); Bournaud et al. (2011). Two items were marked "from memory" in the supplied text (the Hopkins et al. 2016 title and the Volonteri et al. author list). The single most valuable check, if any is run, is eq. 2 and the eps_T range in the 2017 paper, since section 2 of this note rests on it.

### Update (2026-09-21, provenance pass): what was read since

The texts of Angles-Alcazar et al. (2017; arXiv:1603.08007) and Hopkins & Quataert (2011; arXiv:1007.2647) were read (extracted with `pdftotext`). Confirmed in the text: equation 2 of Angles-Alcazar et al. in the form used in section 2 above (the `(M_BH/1e8)^(1/6)`, `(R_0/100 pc)^(-3/2)`, `(1 + f_0/f_gas)^(-1)` dependence, `f_0 = 0.31 f_d^2 (M_d/1e9)^(-1/3)`, and the definition of `f_d`); the normalisation `eps_T` is varied as 0.05, 0.5 and 5 in their runs, with a fiducial 0.5 that they describe as "a factor ~ 10 lower than estimated from simulations without black hole feedback in Hopkins & Quataert (2011)", and it is described as capturing unresolved processes. Hopkins & Quataert state that their analytic accretion-rate predictor was tested against hydrodynamic simulations at a range of galactic scales.

Consequence: "calibrated against hydrodynamical simulations" is supported for the analytic form of Hopkins & Quataert (tested against simulations), not for the normalisation used by Angles-Alcazar et al., which is an adjustable parameter. The manuscript sentence was corrected accordingly. Still not checked: the `eps_m` of about 0.1 and the other items in the list above. The section 2 arithmetic remains an indicative single-nucleus check, not a calibration.