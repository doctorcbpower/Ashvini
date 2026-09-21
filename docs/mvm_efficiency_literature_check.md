# Where do eps_sf = 0.015 and eta_acc = 0.005 come from? A literature check

Purpose: after the C18 tests (crib sheet C18), the black hole-to-stellar mass ratio in the accessible regime is approximately
eta_acc / eps_sf times a feedback factor. Whether that ratio is defensible therefore depends on whether the two efficiencies are
independently motivated. This note records what was checked, what was verified against sources, and what was not. It is a check, not an
experiment, and it changes no result.

Method and limits: the manuscript text and the session notes (`docs/2026_paper_session_code_catalogue.md`, section 6.7) were read
directly. External sources were read through web fetches on 2026-09-20; the fetch tool returns summaries, so quoted equations and table rows
below are "as extracted" and should be checked against the papers. Abstracts only were read for Hobbs et al. (2012) and Hopkins & Quataert
(2011); the full texts were not.

## 1. What the manuscript and the session notes already say

* Manuscript (Section 3, fiducial parameters): "We adopt eps_sf = 0.015, matching the fiducial star-formation efficiency of the underlying
  Ashvini model (Menon, Balu & Power 2026). We set eta_acc = 0.005 as a representative nuclear accretion efficiency ... (Hobbs 2012; Hopkins &
  Quataert 2011)."
* Session notes (catalogue 6.7): eps_sf was moved from 0.02 to 0.015 solely to match the underlying model's published fiducial. eta_acc "has no
  such precedent: it is a new mechanism introduced by this paper, not part of the Ashvini model, so it remains a free parameter justified only
  by the viability margin and the qualitative expectation that nuclear BH accretion is less efficient than star formation given the
  angular-momentum transport bottleneck".
* The "viability margin" is the earlier model's requirement that eta_acc/eps_sf stay below f_BH (memory note "Ashvini viability theorem"). So
  the ratio of the two efficiencies was constrained relative to the target ratio f_BH = 0.5. C18 shows why: in the MVM the
  black-hole-to-stellar ratio of an accessible, supply-limited system relaxes toward roughly eta_acc/eps_sf, so whether the answer is a
  target-sized seed or a small seed depends on whether that ratio (times the feedback factor) lies below or above f_BH.

## 2. What was verified

### eps_sf (star-formation efficiency)

* Menon, Balu & Power (2026), arXiv:2508.08363, "On bursty star formation during cosmological reionization: influence on the metal and dust
  content of low-mass galaxies" (https://arxiv.org/abs/2508.08363). Read from the arXiv HTML version (https://arxiv.org/html/2508.08363), as
  extracted:
  * Table 1 lists eps_sf = 0.015 as the fiducial value (used in their Equations 2 and 5).
  * The star-formation rate is Mdot_star = eps_sf M_g / tau_sf with tau_sf = 0.15 f_sf / H(z), f_sf = 1 "for simplicity", described as a fraction
    of the Hubble time and as "a reasonable upper limit" on the timescale; no observational calibration of eps_sf and no alternative range are
    given; the paper provides minimal direct justification for the value.
* Consequence 1: eps_sf = 0.015 in that paper multiplies a Hubble-time-scale clock, not a local free-fall time. In the MVM the same number
  multiplies M_gas / t_ff(R) at each scale, where t_ff is the free-fall time on the enclosed mass. The value was inherited; its meaning was
  translated. At the galaxy scale the two clocks are within a factor of about 2 (tens of Myr); at the nuclear scale t_ff is of order 1 Myr, so the
  same number implies a specific nuclear star-formation rate about 100 times higher.
* Consequence 2 (minor, irrelevant to the MVM): the code's optional "hubble" clock uses 0.141 t_H(z) (docstring citing Menon & Power 2024,
  Eq. 2), whereas Menon, Balu & Power (2026) as extracted give 0.15/H(z), a difference of about 6 per cent.
* Observational range of the star-formation efficiency per free-fall time (web search, 2026-09-20): a few per cent or less; about 0.01 in the
  theory and observation literature (Krumholz & Tan 2007; McKee & Ostriker 2007; Krumholz et al. 2012), with observed values of order 0.01 to 0.1
  across environments and scales, and some measurements near 0.7 per cent. See for example Utomo et al. (2018),
  https://iopscience.iop.org/article/10.3847/2041-8213/aacf8f and https://arxiv.org/html/1806.11121. As an efficiency per free-fall time,
  0.015 lies within that range. As a nuclear-scale efficiency at 100 to 300 pc, this literature does not say.

### eta_acc (nuclear accretion efficiency)

* Hobbs, Power, Nayakshin & King (2012), MNRAS 421, 3443, "Modelling supermassive black hole growth: towards an improved sub-grid prescription"
  (https://academic.oup.com/mnras/article/421/4/3443/1096548). The abstract argues that Bondi-Hoyle is unreliable as a sub-grid rate and proposes
  an expression that interpolates between the free-fall regime and the Bondi-Hoyle regime. It supports the free-fall form of the supply rate.
  The abstract mentions no efficiency factor and no value of eta_acc.
* Hopkins & Quataert (2011), MNRAS 415, 1027 (https://academic.oup.com/mnras/article/415/2/1027/1033432). The abstract describes an analytic model of
  angular-momentum transport by gravitational torques and "a new estimate of the BH accretion rate given galaxy properties at larger radii, for
  use in galaxy and cosmological simulations and semi-analytic models". The abstract quotes no scaling and no efficiency.
* So the cited sources, as far as their abstracts show, support the form of the nuclear supply (free-fall) and the existence of an
  angular-momentum bottleneck. They do not, in the abstracts, supply the value 0.005 or an order-of-magnitude range for it. The full texts
  were not read, so the stronger statement "they do not contain such a range" is not established.

## 3. Conclusions

1. eps_sf: inherited as a number, with a defensible order of magnitude if read as an efficiency per free-fall time (about 0.01 to 0.1), but
   not derived for the MVM's nuclear scale. Say "adopted from Menon, Balu & Power (2026); the same value is applied per local free-fall time in the MVM".
2. eta_acc: an effective, unresolved feeding parameter. It has no independent calibration in the verified sources. It is a new free parameter
   whose earlier constraint was the requirement that eta_acc/eps_sf lie below f_BH.
3. The two efficiencies are not independently motivated as a pair. The absolute black hole-to-stellar partition in the accessible regime
   is therefore conditional on eta_acc/eps_sf, which is a modelling choice.
4. What the model can still say robustly is how the critical-seed boundary responds to the baryon cycle and to accessibility (claim map H1 to H3,
   R1 to R5). The critical seed being close to the target requires either poor accessibility (fiducial) or eta_acc/eps_sf times the feedback factor
   below f_BH (C18): this dependence should appear in the paper as an explicit condition.

## 4. Not done

* Full texts of Hobbs et al. (2012) and Hopkins & Quataert (2011) were not read; an independent value or range for the fraction of gas at 100 to 300 pc
  that reaches the accretion flow was not located.
* No claim is made that 0.005 is wrong. The check finds that it is a free parameter.

## 5. Update (2026-09-21, provenance pass)

* Full texts were read (extracted with `pdftotext`) for Krumholz & Tan (2007, arXiv:astro-ph/0606277), McKee & Ostriker (2007, arXiv:0707.3514), Hopkins & Quataert (2011, arXiv:1007.2647) and Angles-Alcazar et al. (2017, arXiv:1603.08007). The Hobbs et al. (2012) full text was not.
* Krumholz & Tan: "only ~ 1% of the gas forms stars every free-fall time" in giant molecular clouds, with no evidence of a transition over three orders of magnitude in density. McKee & Ostriker: `eps_ff` "generally low (< ~ 0.01) over a wide range of density tracers", `eps_ff,GMC ~ 0.01`, a theoretical typical value of about 0.02, HCN-traced values of 0.002 to 0.006, and about 0.1 for individual cores. These support "of order a per cent" as used in the manuscript; neither constrains nuclear regions at 100 to 300 pc at z > 7.
* Bibliographic metadata was checked against Crossref: Angles-Alcazar et al. is MNRAS 464, 2840-2853 (2017), DOI 10.1093/mnras/stw2565 (the page number of section 2 above, and the arXiv record's page, were wrong); Krumholz & Tan is ApJ 654, 304-315, DOI 10.1086/509101; McKee & Ostriker is ARA&A 45, 565-687, DOI 10.1146/annurev.astro.45.051806.110602. The existing manuscript bib entry for Hopkins & Quataert (2011) had the wrong arXiv number (1010.1004, a different paper); the correct number is 1007.2647.