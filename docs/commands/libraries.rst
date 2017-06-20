libraries
=========

This section is auto-generated from the help text for the arrow command
``libraries``.


``copy_from_dataset`` command
-----------------------------

**Usage**::

    parsec libraries copy_from_dataset [OPTIONS] LIBRARY_ID DATASET_ID

**Help**

Copy a Galaxy dataset into a library.

**Options**::


      --folder_id TEXT  id of the folder where to place the uploaded files. If not
                        provided, the root folder will be used
      --message TEXT    message for copying action
      -h, --help        Show this message and exit.
    

``create_folder`` command
-------------------------

**Usage**::

    parsec libraries create_folder [OPTIONS] LIBRARY_ID FOLDER_NAME

**Help**

Create a folder in a library.

**Options**::


      --description TEXT     description of the new folder in the data library
      --base_folder_id TEXT  id of the folder where to create the new folder. If not
                             provided, the root folder will be used
      -h, --help             Show this message and exit.
    

``create_library`` command
--------------------------

**Usage**::

    parsec libraries create_library [OPTIONS] NAME

**Help**

Create a data library with the properties defined in the arguments.

**Options**::


      --description TEXT  Optional data library description
      --synopsis TEXT     Optional data library synopsis
      -h, --help          Show this message and exit.
    

``delete_library`` command
--------------------------

**Usage**::

    parsec libraries delete_library [OPTIONS] LIBRARY_ID

**Help**

Delete a data library.

**Options**::


      -h, --help  Show this message and exit.
    

``delete_library_dataset`` command
----------------------------------

**Usage**::

    parsec libraries delete_library_dataset [OPTIONS] LIBRARY_ID DATASET_ID

**Help**

Delete a library dataset in a data library.

**Options**::


      --purged    Indicate that the dataset should be purged (permanently deleted)
      -h, --help  Show this message and exit.
    

``get_folders`` command
-----------------------

**Usage**::

    parsec libraries get_folders [OPTIONS] LIBRARY_ID

**Help**

Get all the folders or filter specific one(s) via the provided ``name`` or ``folder_id`` in data library with id ``library_id``. Provide only one argument: ``name`` or ``folder_id``, but not both.

**Options**::


      --folder_id TEXT  filter for folder by folder id
      --name TEXT       filter for folder by name. For ``name`` specify the full
                        path of the folder starting from the library's root folder,
                        e.g. ``/subfolder/subsubfolder``.
      -h, --help        Show this message and exit.
    

``get_libraries`` command
-------------------------

**Usage**::

    parsec libraries get_libraries [OPTIONS]

**Help**

Get all the libraries or filter for specific one(s) via the provided name or ID. Provide only one argument: ``name`` or ``library_id``, but not both.

**Options**::


      --library_id TEXT  filter for library by library id
      --name TEXT        If ``name`` is set and multiple names match the given name,
                         all the libraries matching the argument will be returned
      --deleted          If set to ``True``, return libraries that have been deleted
      -h, --help         Show this message and exit.
    

``get_library_permissions`` command
-----------------------------------

**Usage**::

    parsec libraries get_library_permissions [OPTIONS] LIBRARY_ID

**Help**

Get the permessions for a library.

**Options**::


      -h, --help  Show this message and exit.
    

``set_library_permissions`` command
-----------------------------------

**Usage**::

    parsec libraries set_library_permissions [OPTIONS] LIBRARY_ID

**Help**

Set the permissions for a library.  Note: it will override all security for this library even if you leave out a permission type.

**Options**::


      --access_in TEXT  list of role ids
      --modify_in TEXT  list of role ids
      --add_in TEXT     list of role ids
      --manage_in TEXT  list of role ids
      -h, --help        Show this message and exit.
    

``show_dataset`` command
------------------------

**Usage**::

    parsec libraries show_dataset [OPTIONS] LIBRARY_ID DATASET_ID

**Help**

Get details about a given library dataset. The required ``library_id`` can be obtained from the datasets's library content details.

**Options**::


      -h, --help  Show this message and exit.
    

``show_folder`` command
-----------------------

**Usage**::

    parsec libraries show_folder [OPTIONS] LIBRARY_ID FOLDER_ID

**Help**

Get details about a given folder. The required ``folder_id`` can be obtained from the folder's library content details.

**Options**::


      -h, --help  Show this message and exit.
    

``show_library`` command
------------------------

**Usage**::

    parsec libraries show_library [OPTIONS] LIBRARY_ID

**Help**

Get information about a library.

**Options**::


      --contents  True if want to get contents of the library (rather than just the
                  library details)
      -h, --help  Show this message and exit.
    

``upload_file_contents`` command
--------------------------------

**Usage**::

    parsec libraries upload_file_contents [OPTIONS] LIBRARY_ID PASTED_CONTENT

**Help**

Upload pasted_content to a data library as a new file.

**Options**::


      --folder_id TEXT  id of the folder where to place the uploaded file. If not
                        provided, the root folder will be used
      --file_type TEXT  Galaxy file format name  [default: auto]
      --dbkey TEXT      Dbkey  [default: ?]
      -h, --help        Show this message and exit.
    

``upload_file_from_local_path`` command
---------------------------------------

**Usage**::

    parsec libraries upload_file_from_local_path [OPTIONS] LIBRARY_ID

**Help**

Read local file contents from file_local_path and upload data to a library.

**Options**::


      --folder_id TEXT  id of the folder where to place the uploaded file. If not
                        provided, the root folder will be used
      --file_type TEXT  Galaxy file format name  [default: auto]
      --dbkey TEXT      Dbkey  [default: ?]
      -h, --help        Show this message and exit.
    

``upload_file_from_server`` command
-----------------------------------

**Usage**::

    parsec libraries upload_file_from_server [OPTIONS] LIBRARY_ID SERVER_DIR

**Help**

Upload all files in the specified subdirectory of the Galaxy library import directory to a library.

**Options**::


      --folder_id TEXT       id of the folder where to place the uploaded files. If
                             not provided, the root folder will be used
      --file_type TEXT       Galaxy file format name  [default: auto]
      --dbkey TEXT           Dbkey  [default: ?]
      --link_data_only TEXT  either 'copy_files' (default) or 'link_to_files'.
                             Setting to 'link_to_files' symlinks instead of copying
                             the files
      --roles TEXT           ???
      -h, --help             Show this message and exit.
    

``upload_file_from_url`` command
--------------------------------

**Usage**::

    parsec libraries upload_file_from_url [OPTIONS] LIBRARY_ID FILE_URL

**Help**

Upload a file to a library from a URL.

**Options**::


      --folder_id TEXT  id of the folder where to place the uploaded file. If not
                        provided, the root folder will be used
      --file_type TEXT  Galaxy file format name  [default: auto]
      --dbkey TEXT      Dbkey  [default: ?]
      -h, --help        Show this message and exit.
    

``upload_from_galaxy_filesystem`` command
-----------------------------------------

**Usage**::

    parsec libraries upload_from_galaxy_filesystem [OPTIONS] LIBRARY_ID

**Help**

Upload a set of files already present on the filesystem of the Galaxy server to a library.

**Options**::


      --folder_id TEXT       id of the folder where to place the uploaded files. If
                             not provided, the root folder will be used
      --file_type TEXT       Galaxy file format name  [default: auto]
      --dbkey TEXT           Dbkey  [default: ?]
      --link_data_only TEXT  either 'copy_files' (default) or 'link_to_files'.
                             Setting to 'link_to_files' symlinks instead of copying
                             the files
      --roles TEXT           ???
      -h, --help             Show this message and exit.
    
