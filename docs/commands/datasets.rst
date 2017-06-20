datasets
========

This section is auto-generated from the help text for the arrow command
``datasets``.


``download_dataset`` command
----------------------------

**Usage**::

    parsec datasets download_dataset [OPTIONS] DATASET_ID

**Help**

Download a dataset to file or in memory.

**Options**::


      --file_path TEXT        If this argument is provided, the dataset will be
                              streamed to disk at that path (should be a directory
                              if use_default_filename=True). If the file_path
                              argument is not provided, the dataset content is
                              loaded into memory and returned by the method (Memory
                              consumption may be heavy as the entire file will be in
                              memory).
      --use_default_filename  If this argument is True, the exported file will be
                              saved as file_path/%s, where %s is the dataset name.
                              If this argument is False, file_path is assumed to
                              contain the full file path including the filename.
                              [default: True]
      --wait_for_completion   If this argument is True, this method call will block
                              until the dataset is ready. If the dataset state
                              becomes invalid, a DatasetStateException will be
                              thrown.
      --maxwait FLOAT         Time (in seconds) to wait for dataset to complete. If
                              the dataset state is not complete within this time, a
                              DatasetTimeoutException will be thrown.  [default:
                              12000]
      -h, --help              Show this message and exit.
    

``show_dataset`` command
------------------------

**Usage**::

    parsec datasets show_dataset [OPTIONS] DATASET_ID

**Help**

Get details about a given dataset. This can be a history or a library dataset.

**Options**::


      --deleted        Whether to return results for a deleted dataset
      --hda_ldda TEXT  Whether to show a history dataset ('hda' - the default) or
                       library dataset ('ldda').  [default: hda]
      -h, --help       Show this message and exit.
    
