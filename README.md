# BioTools Package

## Overview
BioTools is a Python package designed for biological data processing and analysis. It provides functionalities for querying UniProt IDs based on gene names, parsing MITAB files, and retrieving protein information from various APIs.

## Installation

```
pip install git+https://github.com/skewer33/BioTools.git
```


## Usage
After installation, you can import the package and use its functionalities as follows:

```python
from BioTools import gene2uniprotid, get_proteins_info, MITAB_parser, savefig
```

## Modules
- **gene2uniprot**: Functions for querying UniProt IDs based on gene names.
- **MITAB_parser**: Class for parsing MITAB files and extracting relevant information.
- **protein_annotation**: Asynchronous functions to retrieve protein information from UniProt and PDB APIs.
- **wrappers**: Decorator function for saving matplotlib figures.

## Testing
Project contains basic unit tests in the `tests/` folder.

Run tests with:

```bash
python run_tests.py
```

or directly:

```bash
python -m unittest discover -s tests -v
```

### Live (internet) tests
For tests that call real APIs, use:

```bash
python run_live_tests.py
```

Live tests are in `tests/integration/test_gene2uniprot_live.py`:
- `PROBE_GENES`: genes for collecting current API mapping.
- `EXPECTED_UIDS`: expected mapping you can fill after probing.

Recommended flow:
1. Run `python run_live_tests.py`.
2. Check printed `gene -> UniProtID` mapping.
3. Copy values into `EXPECTED_UIDS`.
4. Re-run to get a strict regression check.

## Author
Yakov Mokin - mokinyakov@mail.ru

## License
This project is licensed under the MIT License.
