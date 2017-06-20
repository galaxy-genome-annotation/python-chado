organism
========

This section is auto-generated from the help text for the chakin command
``organism``.


``add_organism`` command
------------------------

**Usage**::

    chakin organism add_organism [OPTIONS] GENUS COMMON ABBR

**Help**

Add a new organism to the Chado database


**Output**


Organism information
   
    
**Options**::


      --species TEXT  The species of the organism
      --comment TEXT  A comment / description
      -h, --help      Show this message and exit.
    

``delete_all_organisms`` command
--------------------------------

**Usage**::

    chakin organism delete_all_organisms [OPTIONS]

**Help**

Delete all organisms


**Output**


None
   
    
**Options**::


      --confirm   Confirm that you really do want to delete ALL of the organisms.
      -h, --help  Show this message and exit.
    

``get_organisms`` command
-------------------------

**Usage**::

    chakin organism get_organisms [OPTIONS]

**Help**

Get all organisms


**Output**


Organisms information
   
    
**Options**::


      --genus TEXT    genus filter
      --species TEXT  species filter
      --common TEXT   common filter
      --abbr TEXT     abbr filter
      --comment TEXT  comment filter
      -h, --help      Show this message and exit.
    
