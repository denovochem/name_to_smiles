# Name manipulation

This package contains utilities that automatically attempt to transform or correct input names if the initial name-to-SMILES conversion fails. At the moment there are two major capabilities:

- Peptide shorthand expansion (e.g. `Boc-Ala-Gly-OMe` → `tert-butoxycarbonyl-l-alanyl-glycine methyl ester`)
- Chemical name correction (primarily aimed at OCR/typo artifacts, with validation via OPSIN)


## Peptide shorthand expansion

Applied when a chemical name looks like peptide shorthand, typically containing amino-acid abbreviations delimited by hyphens, such as:

- `Ala-Gly-Ser`
- `Boc-Asp-Lys(Boc)-OMe`
- `cyclo(Ala-Gly-Ser)` or `cyclo[Ala-Gly-Ser]`

## Chemical name correction

Applied when a chemical name looks like it might have OCR/typo artifacts, with validation via OPSIN.

- Examples of artifacts:
  - `3,4-dihydro- 3-hydroxy-4-oxo-l,2,3-benzotriazine` → `3,4-dihydro- 3-hydroxy-4-oxo-1,2,3-benzotriazine`
  - `ethyl perfiuorobutyrate` → `ethyl perfluorobutyrate`
  - `l-mercapto-Z-thiapropane` → `1-mercapto-2-thiapropane`
