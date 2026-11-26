# Design brief

The dashboard is an evidence ledger for the corrected TMS analyses. Its central visual is an estimate and confidence interval crossing a fixed zero line. Every summary table remains available; small p-values do not determine visibility or order.

Palette: mist `#F4F8FC`, ink `#132B42`, blue `#2767A8`, teal `#147D82`, violet `#7159A6`, amber `#94641B`. White surfaces and faint ink rules provide neutral structure. Color identifies evidence roles, not a binary success/failure judgment.

Type: Segoe UI Variable / Segoe UI / system sans for reading and headings; a local monospace stack for exact numeric values, fields, and provenance. No downloaded fonts. Headings are sentence case, with moderate weight and ample space.

Desktop alignment:

```text
fixed narrow rail | title + corrected-input note
                  | headline + zero-centered interval graphic
                  | complete result ledger, controls above table
                  | descriptive QC image and diagnostic tables
                  | source catalog + verification limits
```

On a phone the rail becomes a short top navigation, filters stack, and the ledger scrolls horizontally within its own frame. The page itself must not overflow.

Principles: show uncertainty before thresholding; distinguish participation from information magnitude; name every correction family; preserve missing values as missing; make original fields and source hashes inspectable; separate treatment comparisons from method associations and QC. Focus indicators, semantic labels, and reduced-motion support are required.

Review against brief: a single generic metric card would hide the multiple families and null findings, so the full ledger is the main interaction. The interval graphic supports the scientific question. A significance-only view is optional and off by default. The design needs no framework, external service, patient records, or new analysis.
