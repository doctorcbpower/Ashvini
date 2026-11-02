# The reservoir models of the 2026 papers: what is kept unchanged, what is superseded

`ashvini/` holds three families of code. The first is the general chain that the current work extends (parameters in `config.py`, components in `components/`, time integrators in `steppers/`). The second and third are models written for the 2026 black-hole papers. They are kept as they were so that the published numbers can be reproduced and compared. They do not use the component structure, and they read their parameters when the module is imported, not from a `Config`.

| Module | Status | Used for |
|---|---|---|
| `reservoir_stock.py` | **Kept unchanged.** The Minimal Viable Model (MVM) of the galaxy and nuclear reservoirs. SHA-256 `31c701dd8f298d4b7bfc2bfb4d74f90fbe77106c817670a1432563b86d255f63`, stored in `tests/golden/mvm.json` and checked by `tests/test_characterisation.py`. | "Critical seed masses for massive black holes in the early Universe" (the production run in `scripts/paper_figures/output/`). |
| `reservoir_stock_premvm.py` | **Kept unchanged; superseded.** The exploratory model before the MVM. | Comparison of the MVM with its predecessor on identical trees only (`scripts/shmr_mvm_vs_general.py`). Not used for a published result. |
| `paper_reservoir.py` | **Partly superseded.** Its halo-growth-rate spline, UV and cold/hot suppression, virial-quantity and tree-interpolation helpers are imported by `reservoir_stock.py` and must not change. Its own integrator (`run_reservoir_paper`, the earlier "Differential Growth" model) is superseded by the MVM. | The helpers: the MVM. The integrator: the record. |
| `paper_reservoir_params.py`, `paper_reservoir_params.yaml` | Kept unchanged with `paper_reservoir.py`. | |
| `black_holes_growth_slimdisk.py`, `spin.py` | Active. The slim-disc growth model is also an option of the general chain (`bh.growth_model: "hobbs_slimdisk"`). | |

## Rules

1. **Do not edit the unchanged files.** The hash test fails if `reservoir_stock.py` changes. If a change is needed, make a new module and say in `docs/SUPERSET_DESIGN.md` why the old one could not be changed.
2. **Do not route new work through them.** Physics for the general chain goes in `components/`; a student starting out should not begin from `reservoir_stock.py`.
3. **The MVM could one day be rebuilt from components.** That is a separate decision (`docs/SUPERSET_DESIGN.md`, "Not done in this plan"); until then it stays as it is.

## Where to read more

* `docs/PRODUCTION_PROVENANCE.md`: what generated the production run, what is inferred, what is not reproducible.
* `docs/mvm_numerical_crib_sheet.md`, `docs/mvm_claim_map.md`: the numbers and claims of the paper and where they come from.
* `scripts/paper_figures/SUPERSEDED.md`: exploratory and earlier-manuscript scripts, kept for the record.
* `MODELS.md`, section "Minimal Viable Model": the MVM's equations.
