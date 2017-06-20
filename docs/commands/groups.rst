groups
======

This section is auto-generated from the help text for the arrow command
``groups``.


``add_group_role`` command
--------------------------

**Usage**::

    parsec groups add_group_role [OPTIONS] GROUP_ID ROLE_ID

**Help**

Add a role to the given group.

**Options**::


      -h, --help  Show this message and exit.
    

``add_group_user`` command
--------------------------

**Usage**::

    parsec groups add_group_user [OPTIONS] GROUP_ID USER_ID

**Help**

Add a user to the given group.

**Options**::


      -h, --help  Show this message and exit.
    

``create_group`` command
------------------------

**Usage**::

    parsec groups create_group [OPTIONS] GROUP_NAME

**Help**

Create a new group.

**Options**::


      --user_ids TEXT  A list of encoded user IDs to add to the new group
      --role_ids TEXT  A list of encoded role IDs to add to the new group
      -h, --help       Show this message and exit.
    

``delete_group_role`` command
-----------------------------

**Usage**::

    parsec groups delete_group_role [OPTIONS] GROUP_ID ROLE_ID

**Help**

Remove a role from the given group.

**Options**::


      -h, --help  Show this message and exit.
    

``delete_group_user`` command
-----------------------------

**Usage**::

    parsec groups delete_group_user [OPTIONS] GROUP_ID USER_ID

**Help**

Remove a user from the given group.

**Options**::


      -h, --help  Show this message and exit.
    

``get_group_roles`` command
---------------------------

**Usage**::

    parsec groups get_group_roles [OPTIONS] GROUP_ID

**Help**

Get the list of roles associated to the given group.

**Options**::


      -h, --help  Show this message and exit.
    

``get_group_users`` command
---------------------------

**Usage**::

    parsec groups get_group_users [OPTIONS] GROUP_ID

**Help**

Get the list of users associated to the given group.

**Options**::


      -h, --help  Show this message and exit.
    

``get_groups`` command
----------------------

**Usage**::

    parsec groups get_groups [OPTIONS]

**Help**

Get all (not deleted) groups.

**Options**::


      -h, --help  Show this message and exit.
    

``show_group`` command
----------------------

**Usage**::

    parsec groups show_group [OPTIONS] GROUP_ID

**Help**

Get details of a given group.

**Options**::


      -h, --help  Show this message and exit.
    

``update_group`` command
------------------------

**Usage**::

    parsec groups update_group [OPTIONS] GROUP_ID

**Help**

Update a group.

**Options**::


      --group_name TEXT  A new name for the group. If None, the group name is not
                         changed.
      --user_ids TEXT    New list of encoded user IDs for the group. It will
                         substitute the previous list of users (with [] if not
                         specified)
      --role_ids TEXT    New list of encoded role IDs for the group. It will
                         substitute the previous list of roles (with [] if not
                         specified)
      -h, --help         Show this message and exit.
    
