Chado Library
=============

.. image:: https://travis-ci.org/galaxy-genome-annotation/python-chado.svg?branch=master
    :target: https://travis-ci.org/galaxy-genome-annotation/python-chado

.. image:: https://readthedocs.org/projects/python-chado/badge/?version=latest
    :target: http://python-chado.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

A Python library for interacting with a Chado database.

Installation
------------

.. code:: bash

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

This will create a chakin config file in ~/.chakin.yml

Examples
--------

.. code:: python

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

Or with the Chakin client:

.. code-block:: shell

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

History
-------

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

Scripts
-------

This library additionally ships with a number of useful command line
scripts in the form of a tool called ``chakin``. The documentation covers that in more detail.

License
-------

Available under the MIT License
