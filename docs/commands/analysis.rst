analysis
========

This section is auto-generated from the help text for the chakin command
``analysis``.


``add_analysis`` command
------------------------

**Usage**::

    chakin analysis add_analysis [OPTIONS] NAME PROGRAM PROGRAMVERSION

**Help**

Create an analysis


**Output**


    Analysis information
    
**Options**::


      --date_executed TEXT  analysis date_executed (yyyy-mm-dd)
      -h, --help            Show this message and exit.
    

``get_analyses`` command
------------------------

**Usage**::

    chakin analysis get_analyses [OPTIONS]

**Help**

Get all or some analyses


**Output**


    Analysis information
    
**Options**::


      --name TEXT            analysis name filter
      --program TEXT         analysis program filter
      --programversion TEXT  analysis programversion filter
      --algorithm TEXT       analysis algorithm filter
      --sourcename TEXT      analysis sourcename filter
      --sourceversion TEXT   analysis sourceversion filter
      --sourceuri TEXT       analysis sourceuri filter
      -h, --help             Show this message and exit.
    
