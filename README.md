# Riemann zero statistics — an experimental research program

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22942475.svg)](https://doi.org/10.5281/zenodo.22942475)

A numerical, pre-registered research program on the local statistics of the
nontrivial zeros of the Riemann zeta function and of Dirichlet L-functions:
how prime (power) waves are transmitted into zero gaps and amplitudes, and how
the zero lattice diffracts.

**Author:** Uğur Sezen. **Status:** four research notes written and internally
audited; not yet peer reviewed. Nothing here claims to bear on a proof of the
Riemann Hypothesis — the results are measured laws, stated with their error bars
and their failed alternatives.

## Research notes (PDF)

| # | Title | File |
|---|---|---|
| 1 | The joint law of zero gaps and local maxima of Hardy's Z-function | [`qm_riemann/arxiv_gap_amplitude.pdf`](qm_riemann/arxiv_gap_amplitude.pdf) |
| 2 | The prime-wave anatomy of the gap–amplitude law | [`qm_riemann/arxiv_prime_wave_anatomy.pdf`](qm_riemann/arxiv_prime_wave_anatomy.pdf) |
| 3 | A response theory for the Riemann zero gas | [`qm_riemann/arxiv_response_theory.pdf`](qm_riemann/arxiv_response_theory.pdf) |
| 4 | The Riemann zero lattice as a warm crystal | [`qm_riemann/arxiv_warm_crystal.pdf`](qm_riemann/arxiv_warm_crystal.pdf) |

LaTeX sources sit next to the PDFs.

## How the work is done

- **Pre-registration.** Measurement tasks are numbered (scripts `qm_riemann/NNN_*.py`,
  tasks 1–198). In the later part of the program the measurement tasks have a
  pencil file (`qm_riemann/KALEM_*.md`, 34 of them) that freezes hypotheses,
  thresholds and death conditions *before* the data are looked at, together
  with a sha256 + timestamp file; the git history of this repository (original
  commit dates preserved) is the public record of that order. Exploratory,
  post-hoc comparisons (e.g. task 196) are labelled as such and are not counted
  as tests.
- **No rescue.** Refuted hypotheses are recorded as refuted; the research log
  keeps the failures next to the successes.
- **Adversarial audits.** Results are audited by independent agents before being
  sealed; derivations are audited before predictions are frozen.
- **Independent replication.** `bagimsiz_dogrulama/` contains a from-scratch
  re-implementation (separate code, raw Odlyzko data) that reproduces the core
  numbers of Notes 1–4; see `bagimsiz_dogrulama/KAPTAN_TEFTISI_23EYL2026.md` for the
  audit of that package.

The internal research log and task files are in Turkish
(`qm_riemann/LITERATUR_TARAMASI_16AGU2026.md` is the running ledger); the notes are
in English.

## AI assistance

The numerical pipelines, audits and drafts were produced with extensive
assistance from Claude (Anthropic); the independent replication package was
written with DeepSeek. The author takes responsibility for all content and
errors.

## Reproducing

- Python 3.9+ with the packages in `requirements.txt`.
- Zero data: A. M. Odlyzko's tables (see `qm_riemann/veri_odlyzko/KAYNAK.txt`);
  the main 2×10⁶-zero set is cached as `qm_riemann/128_odl_zeros6_2e6_zeros.npz`.
  L-function zeros (Note 5 islands) are in `qm_riemann/118b_*_zeros.npz`.
- Scripts are numbered by task (`qm_riemann/NNN_*.py`, `NNN_configs/`). Older
  scripts contain absolute paths of the original machine; `araclar/senkron.sh`
  builds the path bridges, and `araclar/yeniden_kur.sh` regenerates the cached
  intermediate data of tasks 184–187 (~40 min) with a gate against the sealed
  numbers.

## License

Code: MIT (see `LICENSE`). Notes, figures and text: CC BY 4.0.

## Citation

Releases are archived at Zenodo. To cite the program as a whole (all
versions), use

> Uğur Sezen, *Riemann zero statistics — an experimental research program*,
> Zenodo, [doi:10.5281/zenodo.22942475](https://doi.org/10.5281/zenodo.22942475).

Each release also has its own DOI, listed on the Zenodo page (v1.0:
[10.5281/zenodo.22942476](https://doi.org/10.5281/zenodo.22942476)). See also
`CITATION.cff`.
