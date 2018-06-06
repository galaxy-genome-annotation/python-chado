phylogeny
=========

This section is auto-generated from the help text for the chakin command
``phylogeny``.


``add_cvterms`` command
-----------------------

**Usage**::

    chakin phylogeny add_cvterms [OPTIONS]

**Help**

Make sure required cvterms are loaded


**Output**


    created cvterms
    
**Options**::


      -h, --help  Show this message and exit.
    

``gene_families`` command
-------------------------

**Usage**::

    chakin phylogeny gene_families [OPTIONS]

**Help**

Adds an entry in the featureprop table in a chado database for each family a gene belongs to (for use in https://github.com/legumeinfo/lis_context_viewer/).


**Output**


    None
    
**Options**::


      --family_name TEXT  Restrict to families beginning with given prefix
      --nuke              Removes all previous gene families data
      -h, --help          Show this message and exit.
    

``gene_order`` command
----------------------

**Usage**::

    chakin phylogeny gene_order [OPTIONS]

**Help**

Orders all the genes in the database by their order on their respective chromosomes in the gene_order table (for use in https://github.com/legumeinfo/lis_context_viewer/).


**Output**


    None
    
**Options**::


      --nuke      Removes all previous gene ordering data
      -h, --help  Show this message and exit.
    

``load_tree`` command
---------------------

**Usage**::

    chakin phylogeny load_tree [OPTIONS] NEWICK ANALYSIS_ID

**Help**

Load a phylogenetic tree (Newick format) into Chado db


**Output**


    Number of inserted trees
    
**Options**::


      --name TEXT            The name given to the phylotree entry in the database
                             (default=<filename>)
      --xref_db TEXT         The name of the db to link dbxrefs for the trees
                             (default: "null")  [default: null]
      --xref_accession TEXT  The accession to use for dbxrefs for the trees (assumed
                             same as name unless otherwise specified)
      --match_on_name        Match polypeptide features using their name instead of
                             their uniquename
      --prefix TEXT          Comma-separated list of prefix to be removed from
                             identifiers (e.g species prefixes when using loading
                             OrthoFinder output)
      -h, --help             Show this message and exit.
    
