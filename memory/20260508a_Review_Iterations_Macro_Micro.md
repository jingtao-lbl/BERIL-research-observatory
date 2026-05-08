# Review Iterations: Macro-Micro Nutrient Co-cycling

**Date:** May 8, 2026 (early morning, continuation of 20260507h)
**Author:** Jing Tao with Claude
**Type:** Research / Review Iteration
**Project:** `macro_micro_nutrient_cocycling`
**Branch:** `projects/macro_micro_nutrient_cocycling`

---

## Summary

This session addressed reviewer feedback across three `/submit` review cycles (reviews 2→3→4). Each cycle used `--model claude-sonnet-4-6`. The core pattern: each fix round resolved the flagged issues, but the next reviewer pass surfaced new issues or re-flagged partially-addressed ones from a different angle. The `has_N_fixation` column naming issue persisted across all three reviews despite being "fixed" twice.

## Review Cycle Overview

| Review | Commit Before | Commit After | Issues Found | Issues Fixed | New in Next |
|--------|--------------|--------------|:------------:|:------------:|:-----------:|
| 2nd (20260507h) | `d0cb34f` | `12c7d33` | 8 (1C, 2H, 5M/L) | 8 | 6 |
| 3rd | `12c7d33` | `2411450` | 6 (3I, 3N) | 6 | 9 |
| 4th | `2411450` | — | 9 (1C, 2H, 3M, 3L) | pending | — |

C=critical, H=high, I=important, N=nice-to-have, M=medium, L=low.

---

## Complete Issue Tracker

### Issues from 2nd Review → Fixed in `12c7d33`

| # | Severity | Issue | Fix Applied | Re-found? |
|---|----------|-------|-------------|-----------|
| 2.1 | CRITICAL | NB02 permutation uses broad N-fixation (`nifH_pfam` included), produces Z=18.4 instead of REPORT's Z=14.7 | Changed NB02 cell 6 `N_fixation` to `['nifH', 'nifD']`; re-executed; verified Z=14.7 | No — fully resolved |
| 2.2 | HIGH | `has_N_fixation` column in CSV encodes broad definition (5,872) but REPORT describes KO-only (2,746) | Updated src/01 to produce KO-only `has_N_fixation` + broad `has_N_fixation_broad`; added note in NB01 | **YES — 3rd review** said code now mismatches CSV artifact; **4th review** still flags column naming |
| 2.3 | HIGH | NB07 has no figure generation code — only displays pre-existing PNG | Inlined forest plot code from src/08 into NB07 cell 6; re-executed | No — fully resolved |
| 2.4 | MEDIUM | Dead `false_discovery_control` import + `scipy>=1.10` should be `>=1.11` | Removed dead import and `rejected = ...` line; updated requirements.txt to `>=1.11` | No — fully resolved |
| 2.5 | MEDIUM | P. fluorescens listed as nifH-positive but not a known diazotroph | Added †footnote to REPORT table and explanatory sentence | No — fully resolved |
| 2.6 | LOW | 3-species count discrepancy (REPORT: 3,406 vs NB06: 3,409 for soil) | Added clarifying sentence in REPORT Section 6 intro | **Partially** — 4th review (issue 4.7) says the explanation is still not explicit enough |
| 2.7 | LOW | Clarify per-species Fisher test design in Methods | Expanded REPORT Methods bullet with 2×2 table and Stouffer explanation | No — fully resolved (4th review confirms correct) |
| 2.8 | LOW | nifH_pfam row in heatmap Panel A conflicts with "72 pairs" FDR statement | Updated Figure 1 caption to clarify nifH_pfam is reference-only, excluded from 72-pair FDR | No — fully resolved |

### Issues from 3rd Review → Fixed in `2411450`

| # | Severity | Issue | Fix Applied | Re-found? |
|---|----------|-------|-------------|-----------|
| 3.1 | Important | Forest plot narrative says "Strongest in Cyanobacteriota" but Fusobacteriota has higher log-OR (+2.93 vs +1.86) | Added Fusobacteriota as largest effect; qualified Cyanobacteriota as "among phyla with n≥400" | No — fully resolved |
| 3.2 | Important | `has_N_fixation` code–artifact mismatch: src/01 now produces KO-only but CSV has broad | Reverted src/01: broad = `has_N_fixation` (matches CSV), added `has_N_fixation_ko` for KO-only | **YES — 4th review** (issue 4.2) still flags NB01 "N+Metal: 5718" as confusing |
| 3.3 | Important | README Reproduction section omits scripts 06–08 | Appended steps 5–8 to README Reproduction section | No — fully resolved (4th review confirms) |
| 3.4 | Nice-to-have | Add `berdl_notebook_utils` to requirements.txt | Added with comment about Spark-only usage | No — fully resolved (4th review confirms) |
| 3.5 | Nice-to-have | Discuss broad Phenazine phi=0.262 in REPORT | Added paragraph in Results after sensitivity check | No — not re-flagged |
| 3.6 | Nice-to-have | Note env classifier keyword priority as caveat | Added Limitation #8 about keyword priority order | No — not re-flagged |

### Issues from 4th Review (pending)

| # | Severity | Issue | Related to Prior? |
|---|----------|-------|-------------------|
| 4.1 | Critical | NB02 only permutes 4 group pairs; REPORT cites Phz×P Z=4.0 from src/02's 5th pair — no notebook output | **NEW** — not caught in reviews 2 or 3 |
| 4.2 | High | `has_N_fixation` column naming still confusing (NB01 prints "N+Metal: 5718" using broad) | **RECURRING** — same root issue as 2.2 and 3.2; third time flagged |
| 4.3 | High | Fragile PFAM merge loop in src/01 (conditional set_index inside loop) | **NEW** — code style/robustness, not a correctness bug |
| 4.4 | Medium | Permutation p-value floor at N=1000 not noted in REPORT | **NEW** — statistical transparency |
| 4.5 | Medium | "Other" env category (46.8% of species) omitted from REPORT table | **NEW** — completeness |
| 4.6 | Medium | No justification for ≥3 phz gene threshold choice | **NEW** — methods justification |
| 4.7 | Low | 3-species count discrepancy still not explicitly explained | **RECURRING** — same as 2.6; second time flagged despite prior fix |
| 4.8 | Low | NB03 covers 2 analytical steps (core/accessory + phylogenetic) | **NEW** — structural preference |
| 4.9 | Low | Manual BH correction could use statsmodels | **NEW** — code style preference |

---

## Recurring Issues Analysis

### `has_N_fixation` Column Naming (flagged 3 times: 2.2, 3.2, 4.2)

**Root cause:** The CSV artifact was generated by Spark with the broad definition (`nifH + nifD + nifH_pfam` → 5,872 species) in the `has_N_fixation` column. The primary analysis in src/02 recomputes N-fixation from individual gene columns using KO-only (`nifH + nifD` → 2,746 species), so results are always correct. But the derived column in the CSV and its usage in NB01 summary statistics creates a naming ambiguity.

**Fix attempts:**
1. **Review 2 fix (12c7d33):** Changed src/01 to produce KO-only as `has_N_fixation` + broad as `has_N_fixation_broad`. Added NB01 note. → **Problem:** src/01 now produces a different CSV than the committed artifact. Re-running the pipeline would silently change the file.
2. **Review 3 fix (2411450):** Reverted src/01 to match the CSV: broad = `has_N_fixation` (matches artifact), added `has_N_fixation_ko`. Updated NB01 to print both. → **Problem:** NB01 still prints "N+Metal: 5718" using the broad column, which can be confused with the REPORT's primary N×Metal = 2,719.

**Why it keeps recurring:** The reviewer each time correctly identifies a user-facing confusion risk, but each fix only addresses one facet (code vs artifact vs notebook display). The complete fix requires: (a) renaming the CSV column to `has_N_fixation_broad`, (b) adding `has_N_fixation_ko`, (c) updating NB01 to print KO-only as the headline number, and (d) re-generating the CSV from the updated src/01 (requires Spark). Since Spark re-run is impractical in this session, the column naming in the existing CSV cannot be changed — only documented.

### 3-Species Count Discrepancy (flagged 2 times: 2.6, 4.7)

**Root cause:** NB06 cell 2 counts 3,409 soil species with environment assignments, but NB06 cell 4 and the REPORT use n=3,406 for the co-occurrence test (3 species lacked gene family data and were excluded from the join). The REPORT was updated (review 2 fix) with a general clarifying sentence, but the 4th reviewer wanted explicit per-number annotation (e.g., "3 species excluded due to incomplete gene family annotations").

**Why it recurred:** The fix was too vague ("species counts reflect those with both environment assignments and gene family data") rather than explicitly citing the 3-species gap.

---

## Truly New Issues per Review

- **2nd review:** All 8 issues were new (first review cycle after enrichments)
- **3rd review:** 4 genuinely new issues (3.1 Fusobacteriota narrative, 3.3 README steps, 3.4 berdl_notebook_utils, 3.6 env classifier caveat); 1 recurring (3.2 = evolution of 2.2); 1 nice-to-have (3.5 broad Phenazine phi)
- **4th review:** 6 genuinely new issues (4.1 missing Phz×P permutation in NB02, 4.3 PFAM merge loop, 4.4 permutation floor, 4.5 "other" env category, 4.6 phz threshold justification, 4.8 NB03 structure); 2 recurring (4.2 = has_N_fixation naming, 4.7 = count discrepancy); 1 code style (4.9 BH via statsmodels)

**Pattern:** Each review catches ~6 new issues not seen before. The reviewer is non-deterministic — it focuses on different aspects each time. The diminishing-returns point has not been reached after 3 fix cycles, but the severity is trending downward (review 2 had 1 critical correctness bug; review 4's critical is a notebook completeness gap).

---

## Commits This Session

| Hash | Message |
|------|---------|
| `12c7d33` | Address all 8 reviewer issues from second /submit review |
| `2411450` | Address 6 issues from third review |

## Files Changed This Session

| File | Action | Description |
|------|--------|-------------|
| `src/01_extract_gene_families.py` | Modified (×2) | N-fixation column: KO-only→broad+ko→broad+ko (reverted to match CSV) |
| `src/02_cooccurrence_stats.py` | Modified | Removed dead `false_discovery_control` import |
| `notebooks/NB01_gene_family_extraction.ipynb` | Modified (×2) | Updated N-fixation summary display twice |
| `notebooks/NB02_cooccurrence_analysis.ipynb` | Modified | Fixed N_fixation to KO-only; re-executed; Z=14.7 confirmed |
| `notebooks/NB07_forest_plot.ipynb` | Modified | Inlined forest plot code from src/08 |
| `REPORT.md` | Modified (×2) | P. fluorescens footnote, count discrepancy note, Fisher test clarification, nifH_pfam heatmap caption, Fusobacteriota forest plot, broad Phenazine phi note, env classifier Limitation #8 |
| `README.md` | Modified | Added scripts 06-08 to Reproduction section |
| `REVIEW.md` | Regenerated (×2) | Third review (6 issues) → fourth review (9 issues) |
| `requirements.txt` | Modified (×2) | scipy>=1.11; added berdl_notebook_utils |
| `beril.yaml` | Modified | Updated last_session_at |

## Remaining Work

1. Address 4th review issues (or override) — user decision pending
   - Issue 4.1 (critical): Add Phz×P pair to NB02 permutation test
   - Issue 4.2 (high): has_N_fixation naming — at this point, a documentation-only fix may be the pragmatic choice since the CSV can't be regenerated without Spark
   - Issue 4.3 (high): Refactor PFAM merge loop (code robustness, not correctness)
   - Issues 4.4–4.9: Medium/low improvements
2. GitHub authentication (`gh auth login`) still needed before push
3. Final commit + push to `origin` (jingtao-lbl fork)
