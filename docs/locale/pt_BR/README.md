# Portuguese (Brazil) translation notes

These catalogs are written for Brazilian users of Asc-Seurat.

Keep common single-cell analysis terms in English when that is the term most
Brazilian researchers are likely to recognize in class or in the app UI. In
particular, do not translate terms such as:

- `single-cell`
- `scRNA-seq` and `RNA-seq`
- `workflow`
- `quality control` and `QC`
- `clustering`, `cluster`, and `clusters`
- `dataset` and `metadata`
- `UMAP`, `t-SNE`, `PCA`, `embedding`, and `embeddings`
- `trajectory inference`
- `cell-type annotation`
- `doublet`
- `marker`, `markers`, and `cluster markers` when used as analysis terms
- package, method, function, file-format, and app UI names

Use Portuguese for explanatory prose, but preserve code literals, Sphinx roles,
cross-reference targets, package names, function names, UI labels, file names,
and command-line snippets exactly.

The GPL license text is intentionally left in English to avoid creating an
unofficial legal translation.

To refresh the catalogs after editing the English `.rst` files, run these from
the `docs/` directory:

```bash
sphinx-build -b gettext . _build/gettext
sphinx-intl update -p _build/gettext -l pt_BR
sphinx-build -b html -D language=pt_BR . _build/html/pt_BR
```
