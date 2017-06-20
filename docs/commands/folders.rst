folders
=======

This section is auto-generated from the help text for the arrow command
``folders``.


``create_folder`` command
-------------------------

**Usage**::

    parsec folders create_folder [OPTIONS] PARENT_FOLDER_ID NAME

**Help**

Create a folder.

**Options**::


      --description TEXT  folder's description
      -h, --help          Show this message and exit.
    

``delete_folder`` command
-------------------------

**Usage**::

    parsec folders delete_folder [OPTIONS] FOLDER_ID

**Help**

Marks the folder with the given ``id`` as `deleted` (or removes the `deleted` mark if the `undelete` param is True).

**Options**::


      --undelete  If set to True, the folder will be undeleted (i.e. the `deleted`
                  mark will be removed)
      -h, --help  Show this message and exit.
    

``get_permissions`` command
---------------------------

**Usage**::

    parsec folders get_permissions [OPTIONS] FOLDER_ID SCOPE

**Help**

Get the permissions of a folder.

**Options**::


      -h, --help  Show this message and exit.
    

``set_permissions`` command
---------------------------

**Usage**::

    parsec folders set_permissions [OPTIONS] FOLDER_ID

**Help**

Set the permissions of a folder.

**Options**::


      --action TEXT      action to execute, only "set_permissions" is supported.
                         [default: set_permissions]
      --add_ids TEXT     list of role IDs which can add datasets to the folder
      --manage_ids TEXT  list of role IDs which can manage datasets in the folder
      --modify_ids TEXT  list of role IDs which can modify datasets in the folder
      -h, --help         Show this message and exit.
    

``show_folder`` command
-----------------------

**Usage**::

    parsec folders show_folder [OPTIONS] FOLDER_ID

**Help**

Display information about a folder.

**Options**::


      -h, --help  Show this message and exit.
    

``update_folder`` command
-------------------------

**Usage**::

    parsec folders update_folder [OPTIONS] FOLDER_ID NAME

**Help**

Update folder information.

**Options**::


      --description TEXT  folder's description
      -h, --help          Show this message and exit.
    
