# MfAP - Molecular feature Association Pipeline

MfAP is a framework that allows users to calculate and predict DNA- and protein-level consequences of genetic variation in single-exon (intronless) genes.


## Web version

The interactive web application (recommended default) is at https://tim-schubert-mfap.share.connect.posit.cloud/.


## Standalone install

While we encourage use of the online version of MfAP (see above), we provide instructions for local use below:

Requirements: **R ≥ 4.4** and **Python 3.9–3.12** (3.9 recommended; TensorFlow does not support 3.13+).

From the repository root (after `git clone` or unpacking a release archive) run:

```bash
bash mfap.sh
```

Later launches:

```bash
bash mfap.sh --run
```

To reinstall:

```bash
bash mfap.sh --force
```

Note: the repository includes TITER model weights under `titer/model/`, so clones and archives are relatively large.

