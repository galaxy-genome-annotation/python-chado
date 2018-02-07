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
    

``load_tree`` command
---------------------

**Usage**::

    chakin phylogeny load_tree [OPTIONS] NEWICK ANALYSIS_ID

**Help**

Load a phylogenetic tree (Newick format) into Chado db


**Output**


    None
    
**Options**::


      --name TEXT            The name given to the phylotree entry in the database
                             (default=<filename>)
      --xref_db TEXT         The name of the db to link dbxrefs for the trees
                             (default: "null")  [default: null]
      --xref_accession TEXT  The accession to use for dbxrefs for the trees (assumed
                             same as name unless otherwise specified)
      --match_on_name        Match polypeptide features usnig their name instead of
                             their uniquename
      --prefix TEXT          Comma-separated list of prefix to be removed from
                             identifiers (e.g species prefixes when using loading
                             OrthoFinder output)
      -h, --help             Show this message and exit.
    
