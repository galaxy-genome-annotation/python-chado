# Chado Library

[![Test](https://github.com/galaxy-genome-annotation/python-chado/workflows/Test/badge.svg)](https://github.com/galaxy-genome-annotation/python-chado/actions/)
[![Documentation](https://readthedocs.org/projects/python-chado/badge/?version=latest)](http://python-chado.readthedocs.io/en/latest/?badge=latest)

A Python library for interacting with a Chado database.

## Installation

```bash
$ pip install chado

# On first use you'll need to create a config file to connect to the database, just run:

$ chakin init
Welcome to Chado's Chakin! (茶巾)
PGHOST: xxxx
PGDATABASE: xxxx
PGUSER: xxxx
PGPASS:
PGPORT: 5432
PGSCHEMA: public
```

This will create a chakin config file in ~/.chakin.yml

## Examples

```python
from chado import ChadoInstance
ci = ChadoInstance(dbhost="localhost", dbname="chado", dbuser="chado", dbpass="chado", dbschema="public", dbport=5432)

# Create human species
org = ci.organism.add_organism(genus="Homo", species="sapiens", common="Human", abbr="H.sapiens")

# Then display the list of organisms
orgs = ci.organism.get_organisms()

for org in orgs:
    print('{} {}'.format(org.genus, org.species))

# Create an analysis
an = ci.analysis.add_analysis(name="My cool analysis",
                                   program="Something",
                                   programversion="1.0",
                                   algorithm="Google",
                                   sourcename="src",
                                   sourceversion="2.1beta",
                                   sourceuri="http://example.org/",
                                   date_executed="2018-02-03")

# And load some data
ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=orgs[0]['organism_id'])
ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an['analysis_id'], organism_id=orgs[0]['organism_id'])
```

Or with the Chakin client:

```bash
$ my_org=`chakin organism add_organism --species sapiens Homo Human H.sapiens  | jq -r '.organism_id'`

$ chakin organism get_organisms
[
    {
        "organism_id": 1133,
        "genus": "Homo",
        "species": "sapiens",
        "abbreviation": "H.sapiens",
        "common_name": "Human",
        "comment": null
    }
]

# Then load some data
$ my_analysis=`chakin analysis add_analysis \
    "My cool analysis" \
    "Something" \
    "v1.0" \
    "src" | jq -r '.analysis_id'`


$ chakin feature load_fasta \
    --analysis_id $my_analysis \
    --sequence_type contig \
    ./test-data/genome.fa $my_org
```

## History

- 2.3.7
    - Fix loading of expression data when first column header is not empty

- 2.3.6
    - Fix loading of GO terms from GFF

- 2.3.5
    - Fix has_table() calls with recent sqlalchemy versions

- 2.3.4
    - Now requires biopython >=1.78
    - Fixes biopython sequence usage in recent biopython

- 2.3.3
    - Now requires python >= 3.6
    - Better error reporting for blast loader

- 2.3.2
    - Fix interproscan loader only loading the first result of XML v5
    - Fix interproscan loader failing to load IPR by name

- 2.3.1
    - Fix data loading in Tripal database

- 2.3.0
    - Fix non working --re_parent option in fasta loader
    - allow connection using a preformatted url (needed by galaxy tools using pgutil)
    - added loading of Blast and InterProScan data
    - moved `chakin feature load_go` to `chakin load go`
    - fix sequence computing when landmark sequence is available in the db
    - add more options to match features in expression matrix loader (query_type, match_on_name, re_name, skip_missing)

- 2.2.6
    - fix requirement name for psycopg2 (name change for version >=2.8)

- 2.2.5
    - Added support for units in expression loaders
    - Fix error in load_gff when no source is specified

- 2.2.4
    - Fix broken --skip_missing option for load_go

- 2.2.3
    - Throw a warning instead of an exception when a GFF target feature does not exist

- 2.2.2
    - Bug fixes and improvements to the expression module

- 2.2.1
    - Minor release to fix broken package at pypi, no code change

- 2.2.0
    - Added feature.load_go() to load GO annotation (blast2go results)
    - Added feature.get_feature_analyses() to fetch the analyses associated with a feature
    - Added feature.get_feature_cvterms() to fetch the cvterms associated with a feature
    - Added support for biomaterial/expression data (as used by tripal_analysis_expression)
    - New --protein_id_attr option for feature.load_gff()

- 2.1.5
    - bugfix: fix features deletion when deleting an analysis

- 2.1.4
    - bugfix: fix sporadic errors with AnalysisFeature class declaration

- 2.1.3
    - bugfix: make --species a mandatory arg for organism creation
    - bugfix: fix features deletion when deleting an analysis or an organism
    - update chado docker image

- 2.1.2
    - skip whole database schema reflection for simple tasks (analysis and organism management)
    - fix polypeptide creation for genes beginning at position 0
    - fix various small bugs in phylogeny and featureprop loading
    - fix bug in cvterm creation
    - fix crashes in gbk/gff exporters

- 2.1.1
    - newick: remove prefix from node labels too
    - newick: fix errors with named internal nodes

- 2.1
    - auto reflect db schema
    - add phylogeny module
    - load features from fasta
    - load features from gff3
    - load featureprops from tabular file
    - make chakin util commands work when db is offline
    - add unit tests

- 2.0
    - "Chakin" CLI utility
    - Complete package restructure
    - Nearly all functions renamed

## License

Available under the MIT License
