# September 2026 licence re-audit

Prompted by NAR referee 1, who asked whether THEOBROMA perpetuates upstream
licence misassignment. All sources carrying `status: default` in a permissive
tier were re-examined at the data level, distinct from article and code
licensing.

Scope: five sources at `default` with tier_rank <= 1 as of 15 September 2026.
TM-MC was resolved first (documentary basis located, status corrected to
`inferred`, tier unchanged). The remaining four are recorded below.

---

## StreptomeDB — 6,351 attested, 3,023 sole-attestation

- March assignment: CC BY 4.0, rank 1, status `default`
- NAR 3.0 paper (10.1093/nar/gkaa868): CC BY 4.0, covers the article
- Site copyright notice (Pharmaceutical Bioinformatics lab, Freiburg):
  "Any duplication or use of objects such as texts or diagrams, in other
  electronic or printed publications is not permitted without the
  Pharmaceutical Bioinformatics lab prior agreement."
- Citation requested by accession date, time and URL
- Bioregistry: no data licence recorded
- Third-party corroboration: lotusnprod/lotus-processor,
  docs/licenses/streptomedb.md records the same restrictive terms
- OUTCOME: Unspecified, rank 5, status `inferred`
- Note: version ingested is 3.0 (5 October 2024); 4.0 released later

## MycoCentral — 1,603 attested

- March assignment: CC BY 4.0, rank 1, status `default`
- Legal notices page (mycocentral.eu/funding), publisher ANSES:
  "Partial data extraction is allowed after email registration. Use and
  re-distribution of the data, in part, is allowed but not for commercial
  purposes."
- Further conditions: attribution to the MycoCentral Database required in any
  resulting publication or database; data must be reproduced in full without
  alteration or addition; must remain freely downloadable
- This is an explicit statement about the data, not the article
- OUTCOME: CC BY-NC 4.0, rank 2, status `explicit`

## NaturAr — 880 attested

- March assignment: CC BY 4.0, rank 1, status `default`
- Site describes the database as "collaborative and open-source" and "freely
  available online"; no terms of use page, no data licence
- Preprint (10.26434/chemrxiv-2024-56rks) discusses licensing at length but
  only for the source code: GPL chosen for compatibility with OpenBabel, with
  a table of third-party library licences. No statement covers the data.
- The preprint itself is CC BY-NC 4.0
- OUTCOME: Unspecified, rank 5, status `inferred`

## AMDB — 758 attested

- March assignment: CC BY 4.0, rank 1, status `default`
- Site (amdb.online) describes the database as "freely accessible"; no terms
  of use, no data licence, no download conditions
- Article (10.3390/metabo13101088) is CC BY as MDPI standard; covers the paper
- OUTCOME: Unspecified, rank 5, status `inferred`

---

## Summary

Three sources moved from CC BY 4.0 to Unspecified; one moved from CC BY 4.0 to
CC BY-NC 4.0. Each for a different reason: restrictive site terms; explicit
non-commercial data terms; code licensed but data not; access stated but no
terms at all.

Two distinct findings, which should not be conflated:

1. For StreptomeDB, NaturAr and AMDB the March audit correctly recorded that no
   data-level signal was found (`status: default`) but assigned a permissive
   tier, whereas phytochemdb under the same status received a conservative one.
   The defect is an inconsistent default direction, now applied uniformly.

2. For MycoCentral an explicit data-level statement exists and was not located
   in March. This is a missed term, not a policy inconsistency.
