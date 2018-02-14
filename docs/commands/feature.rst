feature
=======

This section is auto-generated from the help text for the chakin command
``feature``.


``delete_features`` command
---------------------------

**Usage**::

    chakin feature delete_features [OPTIONS]

**Help**

Get all or some features


**Output**


    Features information
    
**Options**::


      --organism_id TEXT  organism_id filter
      --analysis_id TEXT  analysis_id filter
      --name TEXT         name filter
      --uniquename TEXT   uniquename filter
      -h, --help          Show this message and exit.
    

``get_features`` command
------------------------

**Usage**::

    chakin feature get_features [OPTIONS]

**Help**

Get all or some features


**Output**


    Features information
    
**Options**::


      --organism_id TEXT  organism_id filter
      --analysis_id TEXT  analysis_id filter
      --name TEXT         name filter
      --uniquename TEXT   uniquename filter
      -h, --help          Show this message and exit.
    

``load_fasta`` command
----------------------

**Usage**::

    chakin feature load_fasta [OPTIONS] FASTA ORGANISM_ID

**Help**

Load features from a fasta file


**Output**


    None
    
**Options**::


      --sequence_type TEXT    Sequence type  [default: contig]
      --analysis_id INTEGER   Analysis ID
      --re_name TEXT          Regular expression to extract the feature name from
                              the fasta sequence id (first capturing group will be
                              used).
      --re_uniquename TEXT    Regular expression to extract the feature name from
                              the fasta sequence id (first capturing group will be
                              used).
      --match_on_name         Match existing features using their name instead of
                              their uniquename
      --update                Update existing feature with new sequence instead of
                              throwing an error
      --db INTEGER            External database to cross reference to.
      --re_db_accession TEXT  Regular expression to extract an external database
                              accession from the fasta sequence id (first capturing
                              group will be used).
      --rel_type TEXT         Relation type to parent feature ('part_of' or
                              'derives_from').
      --re_parent TEXT        Regular expression to extract parent uniquename from
                              the fasta sequence id (first capturing group will be
                              used).
      --parent_type TEXT      Sequence type of the parent feature
      -h, --help              Show this message and exit.
    

``load_gff`` command
--------------------

**Usage**::

    chakin feature load_gff [OPTIONS] GFF ANALYSIS_ID ORGANISM_ID

**Help**

Load features from a gff file


**Output**


    None
    
**Options**::


      -h, --help  Show this message and exit.
    
