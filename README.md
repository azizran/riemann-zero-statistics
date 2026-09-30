# Riemann zero statistics — an experimental research program

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22942475.svg)](https://doi.org/10.5281/zenodo.22942475)

A numerical, pre-registered research program on the local statistics of the
nontrivial zeros of the Riemann zeta function and of Dirichlet L-functions:
how prime (power) waves are transmitted into zero gaps and amplitudes, and how
the zero lattice diffracts.

**Author:** Uğur Sezen. **Status:** six research notes written and internally
audited (Notes 1–4, 7 and 8; Notes 5 and 6 are in preparation); none of the notes
has been peer reviewed yet. Nothing here claims to bear on a proof of the
Riemann Hypothesis — the results are measured laws, stated with their error bars
and their failed alternatives.

## Research notes (PDF)

| # | Title | File |
|---|---|---|
| 1 | The joint law of zero gaps and local maxima of Hardy's Z-function | [`qm_riemann/arxiv_gap_amplitude.pdf`](qm_riemann/arxiv_gap_amplitude.pdf) |
| 2 | The prime-wave anatomy of the gap–amplitude law | [`qm_riemann/arxiv_prime_wave_anatomy.pdf`](qm_riemann/arxiv_prime_wave_anatomy.pdf) |
| 3 | A response theory for the Riemann zero gas | [`qm_riemann/arxiv_response_theory.pdf`](qm_riemann/arxiv_response_theory.pdf) |
| 4 | The Riemann zero lattice as a warm crystal | [`qm_riemann/arxiv_warm_crystal.pdf`](qm_riemann/arxiv_warm_crystal.pdf) |
| 7 | Arithmetic satellites of the Bragg comb of the Riemann zero lattice: residue classes, characters, and a blind test of the Bogomolny–Keating mirror law | [`qm_riemann/arxiv_comb_satellites.pdf`](qm_riemann/arxiv_comb_satellites.pdf) — preprint [doi:10.5281/zenodo.22960606](https://doi.org/10.5281/zenodo.22960606) |
| 8 | The hump between two close zeros: random unitary matrices and the Riemann zeta function | [`qm_riemann/arxiv_small_gap_hump.pdf`](qm_riemann/arxiv_small_gap_hump.pdf), supplementary material [`qm_riemann/arxiv_small_gap_hump_supp.pdf`](qm_riemann/arxiv_small_gap_hump_supp.pdf) — preprint [doi:10.5281/zenodo.23010708](https://doi.org/10.5281/zenodo.23010708) |

LaTeX sources sit next to the PDFs.

## How the work is done

- **Pre-registration.** Measurement tasks are numbered (scripts `qm_riemann/NNN_*.py`,
  tasks 1–201 and 203). In the later part of the program the measurement tasks have a
  pencil file (`qm_riemann/KALEM_*.md`, 38 of them) that freezes hypotheses,
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

Tasks 200-C and 201 also use zeros from the LMFDB database of D. Platt (four files
`zeros_8846000.dat`, `zeros_99146000.dat`, `zeros_997946000.dat`, `zeros_30599546000.dat`,
the first 1.5·10⁶ zeros of each). They are not redistributed here; download them from the
LMFDB data pages into `qm_riemann/veri_lmfdb/` — their sha256 stamps are in
`qm_riemann/200_configs/ONKAYIT_200C.json`.

Task 203 (`qm_riemann/203_configs/`, report `qm_riemann/203_ORAN_KAPPA_RAPOR.md`) evaluates the
ratios conjecture of Conrey, Farmer and Zirnbauer for the variance and third cumulant of
log|ζ(1/2+it)| and log|ζ'(ρ)| at finite height (Note 8, Section 9 and Appendix A; figure script
`qm_riemann/figures/ladder_ratios.py`). The predictions were committed in the author's private
working repository before they were compared with the sealed data of tasks 200-A/C; the commit
hashes quoted in Note 8 refer to that repository. The check against zeros
(`203_configs/sifir_sinamasi_b1k3/`) also needs the LMFDB file `zeros_8846000.dat`.

## License

Code: MIT (see `LICENSE`). Notes, figures and text: CC BY 4.0.

## Citation

Releases are archived at Zenodo. To cite the program as a whole (all
versions), use

> Uğur Sezen, *Riemann zero statistics — an experimental research program*,
> Zenodo, [doi:10.5281/zenodo.22942475](https://doi.org/10.5281/zenodo.22942475).

Each release also has its own DOI, listed on the Zenodo page (v1.0:
[10.5281/zenodo.22942476](https://doi.org/10.5281/zenodo.22942476)). Note 7 has its own preprint record:

> Uğur Sezen, *Arithmetic satellites of the Bragg comb of the Riemann zero lattice: residue
> classes, characters, and a blind test of the Bogomolny–Keating mirror law*, Zenodo (2026),
> [doi:10.5281/zenodo.22960606](https://doi.org/10.5281/zenodo.22960606).

Note 8 has its own preprint record:

> Uğur Sezen, *The hump between two close zeros: an exact small-gap law for random unitary
> matrices and a tilt ladder for the Riemann zeta function*, Zenodo (2026),
> [doi:10.5281/zenodo.23010708](https://doi.org/10.5281/zenodo.23010708).

See also
`CITATION.cff`.
