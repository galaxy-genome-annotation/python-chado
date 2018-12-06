expression
==========

This section is auto-generated from the help text for the chakin command
``expression``.


``add_biomaterial`` command
---------------------------

**Usage**::

    chakin expression add_biomaterial [OPTIONS] BIOMATERIAL_NAME ORGANISM_ID

**Help**

Add a new biomaterial to the database


**Output**


    Biomaterial details
    
**Options**::


      --description TEXT           Description of the biomaterial
      --biomaterial_provider TEXT  Biomaterial provider name
      --biosample_accession TEXT   Biosample accession number
      --sra_accession TEXT         SRA accession number
      --bioproject_accession TEXT  Bioproject accession number
      --attributes TEXT            Custom attributes (In JSON dict form)
      -h, --help                   Show this message and exit.
    

``add_expression`` command
--------------------------

**Usage**::

    chakin expression add_expression [OPTIONS] ORGANISM_ID ANALYSIS_ID

**Help**

Add an expression matrix file to the database


**Output**


    I have no idea
    
**Options**::


      --separator TEXT  Separating character in the matrix file (ex : ','). Default
                        character is tab.  [default:        ]
      -h, --help        Show this message and exit.
    

``delete_all_biomaterials`` command
-----------------------------------

**Usage**::

    chakin expression delete_all_biomaterials [OPTIONS]

**Help**

Delete all biomaterials


**Output**


    None
    
**Options**::


      --confirm   Confirm that you really do want to delete ALL of the biomaterials.
      -h, --help  Show this message and exit.
    

``delete_biomaterial`` command
------------------------------

**Usage**::

    chakin expression delete_biomaterial [OPTIONS]

**Help**




**Output**


    I have no idea
    
**Options**::


      --names TEXT        JSON list of biomaterial names to delete.  [default: []]
      --ids TEXT          JSON list of biomaterial ids to delete.  [default: []]
      --organism_id TEXT  Delete all biomaterial associated with this organism id.
      --analysis_id TEXT  Delete all biomaterial associated with this analysis id.
      -h, --help          Show this message and exit.
    

``delete_biomaterials`` command
-------------------------------

**Usage**::

    chakin expression delete_biomaterials [OPTIONS]

**Help**

Will delete biomaterials based on selector. Only one selector will be used.


**Output**


    Number of deleted biomaterials
    
**Options**::


      --names TEXT        JSON list of biomaterial names to delete.
      --ids TEXT          JSON list of biomaterial ids to delete.
      --organism_id TEXT  Delete all biomaterial associated with this organism id.
      --analysis_id TEXT  Delete all biomaterial associated with this analysis id.
      -h, --help          Show this message and exit.
    

``get_biomaterials`` command
----------------------------

**Usage**::

    chakin expression get_biomaterials [OPTIONS]

**Help**

List biomaterials in the database


**Output**


    List of biomaterials
    
**Options**::


      --provider_id TEXT       Limit query to the selected provider
      --biomaterial_id TEXT    Limit query to the selected biomaterial id
      --organism_id TEXT       Limit query to the selected organism
      --biomaterial_name TEXT  Limit query to the selected biomaterial name
      --analysis_id TEXT       Limit query to the selected analysis_id
      -h, --help               Show this message and exit.
    
