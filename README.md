# Chado Library

A Python library for interacting with a Chado database.

## Scripts

This library additionally ships with a number of useful command line scripts:

Script               | Description
-------------------- | ------------
`create_analysis.py` | Create a new analysis
`create_organism.py` | Create a new organism
`dbshell.py`         | Invoke a database shell session
`export_fa.py`       | Export any sequence records associated with an organism
`export_gbk.py`      | Export a GenBank formatted dataset for an organism
`export_gff3.py`     | Export a GFF3 formatted dataset for an organism
`list_organisms.py`  | List organisms in the database
`purge_organism.py`  | Remove an organism from the database

(This list was generated with `egrep "argparse.*description='(.*)'" -R scripts -o | sed "s/'$//g;s/:.*'/ | /g;"`)

### Optional Dependency: Progress

If the library `tqdm` is installed, a progress bar will be displayed for the
scripts which export data from the database.

## License

Available under the MIT License
