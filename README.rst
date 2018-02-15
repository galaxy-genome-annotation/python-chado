Chado Library
=============

.. image:: https://travis-ci.org/galaxy-genome-annotation/python-chado.svg?branch=master
    :target: https://travis-ci.org/galaxy-genome-annotation/python-chado

.. image:: https://readthedocs.org/projects/python-chado/badge/?version=latest
    :target: http://python-chado.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

A Python library for interacting with a Chado database.

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

Or with the Chakin client:

.. code-block:: shell

    $ chakin organism add_organism --species sapiens Homo Human H.sapiens
    {
        "abbreviation": "H.sapiens",
        "comment": null,
        "common_name": "Human",
        "genus": "Homo",
        "organism_id": 1133,
        "species": "sapiens"
    }
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


History
-------

- 2.1
    - auto reflect db schema
    - add phylogeny module
    - load features from fasta
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
