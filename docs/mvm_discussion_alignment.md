# Discussion 5.1 to 5.3 against the MVM claim map

Scope caveat: this note maps the manuscript's Discussion subsections (5.1 to 5.3 as numbered in the current `.tex`) against the claim map. The author's "5.1 to 5.3" meant three conceptual points on the feeding-efficiency assumptions; that is `docs/mvm_feeding_efficiency_discussion.md`. This note remains valid as a list of old-model statements that conflict with the MVM.

Working note, not manuscript text. The `.tex` is untouched; the wording below is offered for you to accept, change or discard. Item IDs (H, S, R, Q, M) are those of `docs/mvm_claim_map.md`; C-numbers are those of `docs/mvm_numerical_crib_sheet.md`. Status of the numerical campaign: stopped after C18 (committed, `ebf7a44`). No literature round was run for this note. Two things from the literature enter only as stated limits: (i) the efficiency provenance check (`docs/mvm_efficiency_literature_check.md`, where only abstracts of Hobbs et al. 2012 and Hopkins & Quataert 2011 were read), and (ii) the gravitational-torque feeding scaling in section 4, which is taken from your message and has not been checked by me against the source.

The current Discussion still describes the pre-MVM model. Where a sentence rests on a number or a mechanism of that model, it is flagged below.

## 1. The three levels (proposed structure, your wording)

1. **What the model establishes.** Given a nuclear gas supply and an effective BH feeding prescription, the critical initial BH mass required to attain a specified M_BH/M_star,tot is predominantly set by the stellar mass formed and by how far the BH can access the nuclear reservoir. Supported by H1, H2 (M_seed,crit = f_BH M_star,tot / G_BH, G_BH about 1.03 to 1.08), H3, R5.
2. **What the stress tests establish.** The BH-to-stellar ratio reached by light seeds in accessible systems is not an independent characteristic scale: it follows approximately from the competition between BH feeding, star formation and feedback. Supported by S3, S5, R8 (C18). Two qualifications to keep with it: it is finite-time relaxation to z = 5, not an asymptote, and it is absent at 3e13 and at the fiducial accessibility.
3. **What remains unresolved.** The effective BH feeding efficiency is not independently calibrated. Gravitational-torque prescriptions are one physically motivated alternative, but their applicability and calibration in the high-redshift, gas-rich regime are uncertain. Supported by R8 and the provenance note.

## 2. Paragraph by paragraph

### 5.1 What the critical boundary tells us

| Sentence or claim | Status | Reason and suggested direction |
|---|---|---|
| The boundary reflects "the finite time available for Eddington-limited growth ... whether that gas can be converted into black hole mass rapidly enough" | **Conflicts with H2, H3** | At the fiducial the Eddington cap never binds (0% of steps; 0.028 to 0.074 of about 20.8 e-folds used) and M_seed,crit is within 3 to 8% of f_BH M_star,tot. The MVM boundary is set by the stellar mass formed and by access to the gas, not by e-folding time. The e-folding argument applies to the light-seed and feedback-limited regimes (S1), not to the critical seed. This is the largest single change. |
| "Complementary to previous work ... M_seed,crit provides a common quantity" | Compatible | Keep as a conditional statement (the agreed central sentence). |
| Assembly histories give a distribution; nuclear supply and feedback shift the boundary "substantially" | Partly | Assembly: R2 (6 to 10%), Q4 (0.057 to 0.195 dex, defined as a 16-84 width). "Substantially" holds for accessibility through sigma_j (x0.11 to 0.40), not for R_nuc (x0.87 to 1.11), eta_acc (x0.5 to 1.27) or AGN feedback at the critical seed. Say which. |
| Comparison with strict-Eddington versus super-Eddington models | Compatible, needs the H3 qualification | The MVM strict cap is a constraint that does not bind at the critical seed; it binds for light seeds and in early feedback-limited steps (C14, Q3). |
| "Models in which rapid growth drives the ratio towards an attractor can erase ... the seed (Hu2025)" | **Needs care** | C18 finds restoring behaviour in the accessible regime, largely mechanical (competing sinks in one reservoir, d ln R/dt = Mdot_BH/M_BH - Mdot_star/M_star), with a normalisation set by eta_acc/eps_sf and feedback and finite-time relaxation to z = 5. Do not describe the MVM as having an attractor; do say that seed memory is strongly reduced in the accessible regime at 3e10 and 3e11 and not at 3e13 or at the fiducial. I have not read Hu et al. (2025) and make no statement about their model. |
| "Interpreted as the boundary for the adopted Eddington-limited gas-accretion pathway" | Keep, with one addition | Add that the pathway is defined by the effective feeding prescription, whose efficiency is uncalibrated (level 3). |

### 5.2 Why the critical mass is not universal

| Sentence or claim | Status | Reason and suggested direction |
|---|---|---|
| Tree-to-tree spread "0.12 to 0.31 dex" | **Old-model number** | Q4: 0.057 to 0.195 dex for the MVM, as a 16-84 width; not comparable. Crib sheet section 3 says the old numbers must not appear. |
| Nuclear supply and feedback "move the boundary by factors of tens to more than a hundred" | **Not supported by the MVM** | R5, R6: sigma_j 1.5 gives x0.11 to 0.14 (about a factor 9); R_nuc x0.87 to 1.11; eta_acc x0.1 / x10 gives 1.17 to 1.27 / 0.52 to 0.75. The MVM does not reproduce a factor of tens to a hundred from these parameters. |
| "These parameters ... represent physical processes whose efficiency is not yet independently calibrated" | **Supported and strengthened** | R8 and the provenance note: eps_sf inherited (a Hubble-scale clock in the source, applied per local free-fall time here); eta_acc has no independent calibration in the sources checked; the pair is not independently motivated. |
| Radiative efficiency gives "a factor of about 3 to 4" in M_seed,crit | **Not supported by the MVM** | R6: x1.04 to 1.05 (0.057) and x0.87 to 0.91 (0.32), weak because the cap does not bind at the critical seed. |
| Sub-100 pc transport is unresolved; spin and radiative efficiency are prescribed | Keep | This is the natural place for the torque-feeding limitation (section 4 below). |

### 5.3 Implications for high-z overmassive black holes

| Sentence or claim | Status | Reason and suggested direction |
|---|---|---|
| The boundary as a mapping from host properties to required seed | Compatible, conditional | Level 1 holds; the mapping to an observed population additionally needs the feeding efficiency (level 3). |
| "A seed below the boundary does not rule out a channel" | Keep, and sharpen | S1 to S3: light seeds reach 4e-7 to 4e-2 of M_star,tot at the fiducial (H4), and 0.07 to 0.12 at sigma_j = 1.5; conclusion depends on accessibility. |
| Population-level comparison is more informative | Keep | Compatible with Q4 and R2. |
| Accretion-state remarks (GN-z11; Juodzbalis et al.) | Not tested by the MVM | Left as is; nothing in C1 to C18 bears on them. |

Adjacent, outside 5.1 to 5.3 but affected: Conclusions items 2 and 3 ("existence of the boundary is robust ... two decades in nuclear-supply efficiency"; "factors of tens to about 100") and the seed-channel item are old-model statements. The agreed position is not to call the boundary simply "robust": it is insensitive to seed epoch and to stochastic versus smooth assembly, and its high-mass shape depends on the imposed M_hot (R1 to R3).

## 3. What the stress tests add to the discussion, in the model's own terms

* The MVM erases seed memory in accessible systems for a reason that does not require the constant-efficiency construction: the nuclear supply rate is eta_acc M_nuc / t_ff, with t_ff set by the enclosed mass, so M_BH enters only weakly while M_BH is small compared with the nuclear gas and stars. The level of the ratio, not the loss of seed memory, is what depends on eta_acc/eps_sf (C18, R8).
* Accessibility matters for whether light seeds enter rapid growth (H4, R5, S3). This is the one physical dependence with a large lever arm on M_seed,crit in the MVM.
* The 0.1 to 0.2 level is a consequence of the adopted efficiencies and feedback, not a characteristic scale (S3; not-claimed items 13 to 15).

## 4. Where gravitational-torque feeding belongs

As stated in your message (Angles-Alcazar and Hopkins-Quataert framework; the scaling Mdot_BH proportional to f_d^{5/2} M_BH^{1/6} M_d R_0^{-3/2}; not verified by me against the sources):

* It is a different prescription with resolved-variable dependence (disc gas fraction, BH mass, nuclear mass distribution) that the MVM does not have. Its efficiency parameter is not eta_acc, so no sweep in the MVM should be motivated by its literature range: that would conflate two prescriptions.
* It offers a physical route to replacing the constant eta_acc, and its weak M_BH dependence is consistent with seed memory being erased when supply is torque-regulated. In the MVM this is an inference by analogy, not a tested result.
* Appropriate placement: the sub-100 pc transport paragraph in 5.2 (as the source of the unresolved feeding efficiency) and, if wanted, the future-work sentence. The future question is whether replacing the phenomenological constant with a physically motivated nuclear transport law changes the critical-seed result qualitatively. It is not a claim of this paper.

## 5. Not established, so not to be said

* That eta_acc = 0.005 is calibrated by Hobbs et al. (2012) or Hopkins & Quataert (2011); only their abstracts were read, and they do not give the value.
* That the torque model would reproduce, or reverse, any MVM result.
* That the MVM has an attractor, a ceiling or a characteristic M_BH/M_star.
* That Eddington-limited e-folding sets the critical seed at the fiducial.
