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
    

``get_feature_analyses`` command
--------------------------------

**Usage**::

    chakin feature get_feature_analyses [OPTIONS] FEATURE_ID

**Help**

Get analyses associated with a feature


**Output**


    Feature analyses
    
**Options**::


      -h, --help  Show this message and exit.
    

``get_feature_cvterms`` command
-------------------------------

**Usage**::

    chakin feature get_feature_cvterms [OPTIONS] FEATURE_ID

**Help**

Get cvterms associated with a feature


**Output**


    Feature cvterms
    
**Options**::


      -h, --help  Show this message and exit.
    

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


    Number of inserted sequences
    
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
    

``load_featureprops`` command
-----------------------------

**Usage**::

    chakin feature load_featureprops [OPTIONS] TAB_FILE ANALYSIS_ID

**Help**

Load feature properties from a tabular file (Column1: feature name or uniquename, Column2: property value)


**Output**


    Number of inserted featureprop
    
**Options**::


      --feature_type TEXT  Type of the target features in sequence ontology (will
                           speed up loading if specified)
      --match_on_name      Match features using their name instead of their
                           uniquename
      -h, --help           Show this message and exit.
    

``load_gff`` command
--------------------

**Usage**::

    chakin feature load_gff [OPTIONS] GFF ANALYSIS_ID ORGANISM_ID

**Help**

Load features from a gff file


**Output**


    None
    
**Options**::


      --landmark_type TEXT       Type of the landmarks (will speed up loading if
                                 provided, e.g. contig, should be a term of the
                                 Sequence ontology)
      --re_protein TEXT          Replacement string for the protein name using
                                 capturing groups defined by --re_protein_capture
      --re_protein_capture TEXT  Regular expression to capture groups in mRNA name
                                 to use in --re_protein (e.g. "^(.*?)-R([A-Z]+)$",
                                 default="^(.*?)$")  [default: ^(.*?)$]
      --fasta TEXT               Path to a Fasta containing sequences for some
                                 features. When creating a feature, if its sequence
                                 is in this fasta file it will be loaded. Otherwise
                                 for mRNA and polypeptides it will be computed from
                                 the genome sequence (if available), otherwise it
                                 will be left empty.
      --no_seq_compute           Disable the computation of mRNA and polypeptides
                                 sequences based on genome sequence and positions.
      --quiet                    Hide progress information
      --add_only                 Use this flag if you're not updating existing
                                 features, but just adding new features to the
                                 selected analysis and organism. It will speedup
                                 loading, and reduce memory usage, but might produce
                                 errors in case of already existing feature.
      --protein_id_attr TEXT     Attribute containing the protein uniquename. It is
                                 searched at the mRNA level, and if not found at CDS
                                 level.
      -h, --help                 Show this message and exit.
    

``load_go`` command
-------------------

**Usage**::

    chakin feature load_go [OPTIONS] INPUT ORGANISM_ID ANALYSIS_ID

**Help**

Load GO annotation from a tabular file


**Output**


    Number of inserted GO terms
    
**Options**::


      --query_type TEXT      The feature type (e.g. 'gene', 'mRNA', 'polypeptide',
                             'contig') of the query. It must be a valid Sequence
                             Ontology term.  [default: polypeptide]
      --match_on_name        Match features using their name instead of their
                             uniquename
      --name_column INTEGER  Column containing the feature identifiers (2, 3, 10 or
                             11; default=2).  [default: 2]
      --go_column INTEGER    Column containing the GO id (default=5).  [default: 5]
      --re_name TEXT         Regular expression to extract the feature name from the
                             input file (first capturing group will be used).
      --skip_missing         Skip lines with unknown features or GO id instead of
                             aborting everything.
      -h, --help             Show this message and exit.
    
