load
====

This section is auto-generated from the help text for the chakin command
``load``.


``blast`` command
-----------------

**Usage**::

    chakin load blast [OPTIONS] ANALYSIS_ID ORGANISM_ID INPUT

**Help**

Load a blast analysis, in the same way as does the tripal_analysis_blast module


**Output**


    Number of processed hits
    
**Options**::


      --blastdb TEXT     Name of the database blasted against (must be in the Chado
                         db table)
      --blastdb_id TEXT  ID of the database blasted against (must be in the Chado db
                         table)
      --re_name TEXT     Regular expression to extract the feature name from the
                         input file (first capturing group will be used).
      --query_type TEXT  The feature type (e.g. 'gene', 'mRNA', 'polypeptide',
                         'contig') of the query. It must be a valid Sequence
                         Ontology term.  [default: polypeptide]
      --match_on_name    Match features using their name instead of their uniquename
      --skip_missing     Skip lines with unknown features or GO id instead of
                         aborting everything.
      -h, --help         Show this message and exit.
    

``go`` command
--------------

**Usage**::

    chakin load go [OPTIONS] INPUT ORGANISM_ID ANALYSIS_ID

**Help**

Load GO annotation from a tabular file, in the same way as does the tripal_analysis_go module


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
    

``interpro`` command
--------------------

**Usage**::

    chakin load interpro [OPTIONS] ANALYSIS_ID ORGANISM_ID INPUT

**Help**

Load an InterProScan analysis, in the same way as does the tripal_analysis_intepro module


**Output**


    Number of processed hits
    
**Options**::


      --parse_go         Load GO annotation to the database
      --re_name TEXT     Regular expression to extract the feature name from the
                         input file (first capturing group will be used).
      --query_type TEXT  The feature type (e.g. 'gene', 'mRNA', 'polypeptide',
                         'contig') of the query. It must be a valid Sequence
                         Ontology term.  [default: polypeptide]
      --match_on_name    Match features using their name instead of their uniquename
      --skip_missing     Skip lines with unknown features or GO id instead of
                         aborting everything.
      -h, --help         Show this message and exit.
    
