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

## Author
Mokin Yakov - mokinyakov@mail.ru

## License
This project is licensed under the MIT License.