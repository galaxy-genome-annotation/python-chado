genomes
=======

This section is auto-generated from the help text for the arrow command
``genomes``.


``get_genomes`` command
-----------------------

**Usage**::

    parsec genomes get_genomes [OPTIONS]

**Help**

Returns a list of installed genomes

**Options**::


      -h, --help  Show this message and exit.
    

``install_genome`` command
--------------------------

**Usage**::

    parsec genomes install_genome [OPTIONS]

**Help**

Download and/or index a genome.

**Options**::


      --func TEXT           Allowed values: 'download', Download and index; 'index',
                            Index only  [default: download]
      --source TEXT         Data source for this build. Can be: UCSC, Ensembl, NCBI,
                            URL
      --dbkey TEXT          DB key of the build to download, ignored unless 'UCSC'
                            is specified as the source
      --ncbi_name TEXT      NCBI's genome identifier, ignored unless NCBI is
                            specified as the source
      --ensembl_dbkey TEXT  Ensembl's genome identifier, ignored unless Ensembl is
                            specified as the source
      --url_dbkey TEXT      DB key to use for this build, ignored unless URL is
                            specified as the source
      --indexers TEXT       POST array of indexers to run after downloading
                            (indexers[] = first, indexers[] = second, ...)
      -h, --help            Show this message and exit.
    

``show_genome`` command
-----------------------

**Usage**::

    parsec genomes show_genome [OPTIONS] ID

**Help**

Returns information about build <id>

**Options**::


      --num TEXT    num
      --chrom TEXT  chrom
      --low TEXT    low
      --high TEXT   high
      -h, --help    Show this message and exit.
    
