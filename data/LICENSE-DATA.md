# Licences of the data and model files in this directory

Glass's own source code is under the GNU Affero General Public License v3 or later
(`../LICENSE`). The files here are data, and three different sets of terms apply.

## `hose_model.pkl.gz` — NMR shift model

Built from **nmrshiftdb2** (https://nmrshiftdb.nmr.uni-koeln.de), an open database of
experimental, atom-assigned NMR spectra, by Stefan Kuhn, Christoph Steinbeck and contributors.
The model is a *produced work* derived from that database, so the **nmrshiftdb2 Database
Licence** applies to it as well as to the database. That licence is derived from the Open
Database License (ODbL) 1.0 with additions, and it requires:

- **Attribution** to nmrshiftdb2 wherever the model or its predictions are used publicly.
- **Share-alike**: a derived database or produced work must be offered under the same licence.
- **Section 4.5**: software that cannot serve its purpose without the database must itself carry
  an OSI-approved licence. Glass satisfies this with the AGPL-3.0-or-later.

Anyone hosting Glass publicly must keep the attribution visible and offer this model file for
download. The raw database (`nmrshiftdb2withsignals.sd`, 151 MB) is **not** redistributed here;
download it from the project above and rebuild with `python -m core.hose build`.

Cite: Steinbeck, C.; Krause, S.; Kuhn, S. *J. Chem. Inf. Comput. Sci.* **2003**, 43, 1733–1739.

## `periodic_table.json` — element data

From the **PubChem Periodic Table** (NCBI/NLM), https://pubchem.ncbi.nlm.nih.gov/periodic-table/.
US Government work, public domain. No restrictions; attribution kept as a courtesy. PubChem's
usage policy caps automated requests at 5 per second.

## `benchmarks.json` — our own comparison table

Measured in Glass and compared with published experimental frequencies. Part of this project and
covered by the project licence (`../LICENSE`); the underlying experimental values are facts from
the literature cited in `../THIRD_PARTY.md`.
