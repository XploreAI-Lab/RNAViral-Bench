# Licensing

This repository contains separately licensed components.

Unless otherwise indicated, the benchmarking notebooks, execution
wrappers, analysis scripts, and visualization code are licensed under
the GNU General Public License v3.0 only (GPL-3.0-only).

Parts of the benchmarking framework were adapted from
`sinc-lab/lncRNA-folding`, which is distributed under GPL v3.

`methods/ct2dot.py` is separately licensed under the GNU General
Public License v2.0 only (GPL-2.0-only). It is a standalone
command-line Python port/adaptation of functionality from RNAstructure.

The benchmark invokes `ct2dot.py` only as a separate subprocess and
communicates with it through input/output files.

Third-party RNA secondary structure prediction software and pretrained
models are not redistributed in this repository. They remain subject
to their respective upstream licenses.

Complete license texts are available in the `LICENSES/` directory.