# Third-party components and credits

Glass is built on the following open-source software and public data. Versions are
those present when this file was generated; see `requirements.txt` for constraints.

| Component | Role in Glass | Version | License | Source |
|---|---|---|---|---|
| OPSIN (Open Parser for Systematic IUPAC Nomenclature), Daniel Lowe et al., University of Cambridge | Parses IUPAC/systematic names to structures; provides per-atom locants used for the numbering overlay | 2.9.0 (jar bundled by py2opsin) | MIT | https://github.com/dan2097/opsin — Lowe, D. M.; Corbett, P. T.; Murray-Rust, P.; Glen, R. C. *J. Chem. Inf. Model.* **2011**, 51, 739–753 |
| py2opsin | Ships the OPSIN jar; we only use its file location | 1.2.0 | MIT | https://github.com/JacksonBurns/py2opsin |
| RDKit | Cheminformatics toolkit: structure handling, canonical SMILES/InChI, 2D and 3D coordinates, descriptors | 2026.3.6 | BSD-3-Clause | https://www.rdkit.org |
| InChI Trust library (via RDKit) | Standard InChI / InChIKey used for every identity check | bundled | IUPAC/InChI Trust Licence | https://www.inchi-trust.org |
| thermo (Caleb Bell) | Physical-property estimation | 0.6.1 | MIT | https://github.com/CalebBell/thermo |
| FastAPI | HTTP API | 0.141.1 | MIT | https://fastapi.tiangolo.com |
| uvicorn | ASGI server | 0.52.4 | BSD-3-Clause | https://www.uvicorn.org |
| requests | PubChem HTTP client | 2.34.2 | Apache-2.0 | https://requests.readthedocs.io |
| NumPy | Hessian / normal-mode linear algebra | 2.5.2 | BSD-3-Clause | https://numpy.org |
| SciPy | Special functions for hydrogenic orbital cubes (`core/atom.py`) and geometry optimisation (`core/vibrations.py`) | 1.18.1 | BSD-3-Clause | https://scipy.org |
| 3Dmol.js | Browser 3D viewer and vibrational-mode animation (loaded from cdnjs at runtime, not redistributed) | 2.0.4 | BSD-3-Clause | https://3dmol.csb.pitt.edu — Rego, N.; Koes, D. *Bioinformatics* **2015**, 31, 1322–1324 |
| Ketcher (EPAM Systems) | In-browser structure editor for 'Edit structure' (standalone build with Indigo WebAssembly, vendored unmodified in `static/vendor/ketcher/`, with its Apache-2.0 `LICENSE` and `NOTICE` included there) | 3.17.0 | Apache-2.0 | https://github.com/epam/ketcher |
| tblite (Sebastian Ehlert et al.) | Semiempirical electronic-structure calculations | 0.7.0 | LGPL-3.0 | https://tblite.readthedocs.io — Bannwarth, Ehlert, Grimme *J. Chem. Theory Comput.* **2019**, 15, 1652 |
| nmrshiftdb2 (Stefan Kuhn, Christoph Steinbeck et al.) | Experimental reference data (model component; raw dataset not redistributed) | 2024 export | nmrshiftdb2 Database Licence, derived from ODbL 1.0: attribution, share-alike on derived databases and produced works, and section 4.5 requiring dependent software to be OSI-licensed. Raw database not redistributed here; the model trained from it is, under the same terms (`data/LICENSE-DATA.md`) | Steinbeck, Krause, Kuhn *J. Chem. Inf. Comput. Sci.* **2003**, 43, 1733 |
| SMILES2IUPAC-canonical-base (Knowledgator Engineering) | MT5 seq2seq model generating candidate IUPAC names for structures absent from PubChem; every candidate is gated by an OPSIN InChIKey round-trip | 2024 | Apache-2.0 | https://huggingface.co/knowledgator/SMILES2IUPAC-canonical-base |
| chemical-converters | Loader/inference wrapper for the SMILES→IUPAC model | 0.1.2 | Apache-2.0 | https://pypi.org/project/chemical-converters |
| PyTorch (CPU build) | Runs the SMILES→IUPAC model in the name-generation sidecar | 2.14.1+cpu | BSD-3-Clause | https://pytorch.org |
| Hugging Face transformers | MT5 model and tokenizer implementation for the name generator | 4.48.2 | Apache-2.0 | https://github.com/huggingface/transformers |
| tokenizers, safetensors (Hugging Face) | Tokenisation and weight loading for the name generator | 0.21.4 / 0.8.0 | Apache-2.0 | https://github.com/huggingface |
| PySCF + pyscf-properties (Sun et al.) | Quantum-chemistry calculations (orbitals and NMR shieldings) | 2.14.0 / 0.1.0 | Apache-2.0 | Sun et al. *WIREs Comput. Mol. Sci.* **2018**, 8, e1340 |
| Java runtime (OpenJDK) | Required to run OPSIN; installed from the OS packages, not redistributed | 21 | GPLv2 + Classpath Exception | https://openjdk.org |
| PubChem (NCBI / NLM) | Optional online lookup: database names, IUPAC name, synonyms, CAS RN, experimental melting/boiling points, density, solubility; periodic-table element data (`data/periodic_table.json`, vendored) | PUG REST / PUG View | Public-domain US Government data; usage policy: ≤5 requests/s | https://pubchem.ncbi.nlm.nih.gov — Kim, S. et al. *Nucleic Acids Res.* **2023**, 51, D1373–D1380 |

## Methods implemented in this codebase (not third-party code, but published methods)

| Method | Used for | Reference |
|---|---|---|
| Joback–Reid group contributions | Tm, Tb, critical properties | Joback, K. G.; Reid, R. C. *Chem. Eng. Commun.* **1987**, 57, 233–243 |
| Girolami density rule | Liquid/solid density | Girolami, G. S. *J. Chem. Educ.* **1994**, 71, 962–964 |
| Wildman–Crippen (via RDKit) | log P, molar refractivity | Wildman, S. A.; Crippen, G. M. *J. Chem. Inf. Comput. Sci.* **1999**, 39, 868–873 |
| ¹H shift additivity | Shifts of sp³, aromatic and vinylic protons | Pretsch, E.; Bühlmann, P.; Badertscher, M. *Structure Determination of Organic Compounds*, 4th ed., Springer 2009; Pascual, C.; Meier, J.; Simon, W. *Helv. Chim. Acta* **1966**, 49, 164 |
| ¹³C shift additivity | Carbon shifts | Grant, D. M.; Paul, E. G. *J. Am. Chem. Soc.* **1964**, 86, 2984–2990; Pretsch et al. (above) |
| Group-frequency IR table | Simulated IR bands | Socrates, G. *Infrared and Raman Characteristic Group Frequencies*, 3rd ed., Wiley 2001 |
| Woodward–Fieser rules | UV–Vis λmax of dienes, enones, aryl carbonyls | Woodward, R. B. *J. Am. Chem. Soc.* **1941**, 63, 1123; Fieser, L. F.; Fieser, M. *Steroids*, 1959; Scott, A. I. *Interpretation of the Ultraviolet Spectra of Natural Products*, Pergamon 1964 |
| MMFF94s normal-mode analysis | Vibrational frequencies, mode vectors, approximate IR intensities from fixed MMFF charges | Halgren, T. A. *J. Comput. Chem.* **1996**, 17, 490–519 (MMFF94); Wilson, E. B.; Decius, J. C.; Cross, P. C. *Molecular Vibrations*, McGraw-Hill 1955 |
| Displaced-harmonic-oscillator Franck–Condon model | Photoluminescence panel and interactive Franck–Condon diagram (Poisson vibronic intensities e^−ˢSⁿ/n!, Huang–Rhys factor from structural rigidity, mirror-image emission, Stokes shift) | Franck, J. *Trans. Faraday Soc.* **1926**, 21, 536; Condon, E. U. *Phys. Rev.* **1928**, 32, 858; Huang, K.; Rhys, A. *Proc. R. Soc. A* **1950**, 204, 406; Lakowicz, J. R. *Principles of Fluorescence Spectroscopy*, 3rd ed., Springer 2006 |
| ETKDGv3 conformer embedding | 3D embedding before MMFF minimisation | Riniker, S.; Landrum, G. A. *J. Chem. Inf. Model.* **2015**, 55, 2562–2574 |
| ACS 1996 depiction style | 2D structure conventions | *The ACS Style Guide*, 3rd ed., Oxford Univ. Press 2006 |

## Licence of this project

Glass is free software: Copyright (C) 2026 Gregory Pitch, licensed under the **GNU Affero
General Public License, version 3 or later** (`LICENSE`). You may use, study, modify, and
redistribute it, including commercially. If you distribute it, or run a modified version as a
network service, you must offer the corresponding source under the same licence.

This choice is required as well as intended: the NMR shift model ships with data derived from
nmrshiftdb2, whose licence (section 4.5) requires software that depends on it to carry an
OSI-approved licence.

The data and model files have their own terms; see `data/LICENSE-DATA.md`. The third-party
components credited above remain under their own licences.
