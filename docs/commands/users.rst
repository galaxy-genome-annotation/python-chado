users
=====

This section is auto-generated from the help text for the arrow command
``users``.


``create_local_user`` command
-----------------------------

**Usage**::

    parsec users create_local_user [OPTIONS] USERNAME USER_EMAIL PASSWORD

**Help**

Create a new Galaxy local user.

**Options**::


      -h, --help  Show this message and exit.
    

``create_remote_user`` command
------------------------------

**Usage**::

    parsec users create_remote_user [OPTIONS] USER_EMAIL

**Help**

Create a new Galaxy remote user.

**Options**::


      -h, --help  Show this message and exit.
    

``create_user`` command
-----------------------

**Usage**::

    parsec users [OPTIONS] COMMAND [ARGS]...

**Help**

Deprecated method.


``create_user_apikey`` command
------------------------------

**Usage**::

    parsec users create_user_apikey [OPTIONS] USER_ID

**Help**

Create a new API key for a given user.

**Options**::


      -h, --help  Show this message and exit.
    

``delete_user`` command
-----------------------

**Usage**::

    parsec users delete_user [OPTIONS] USER_ID

**Help**

Delete a user.

**Options**::


      --purge     if ``True``, also purge (permanently delete) the history
      -h, --help  Show this message and exit.
    

``get_current_user`` command
----------------------------

**Usage**::

    parsec users get_current_user [OPTIONS]

**Help**

Display information about the user associated with this Galaxy connection.

**Options**::


      -h, --help  Show this message and exit.
    

``get_user_apikey`` command
---------------------------

**Usage**::

    parsec users get_user_apikey [OPTIONS] USER_ID

**Help**

Get the current API key for a given user. This functionality is available since Galaxy ``release_17.01``.

**Options**::


      -h, --help  Show this message and exit.
    

``get_users`` command
---------------------

**Usage**::

    parsec users get_users [OPTIONS]

**Help**

Get a list of all registered users. If ``deleted`` is set to ``True``, get a list of deleted users.

**Options**::


      --deleted TEXT
      --f_email TEXT  filter for user emails. The filter will be active for non-
                      admin users only if the Galaxy instance has the
                      ``expose_user_email`` option set to ``True`` in the
                      ``config/galaxy.ini`` configuration file. This parameter is
                      silently ignored for non-admin users in Galaxy
                      ``release_15.01`` and earlier.
      --f_name TEXT   filter for user names. The filter will be active for non-admin
                      users only if the Galaxy instance has the ``expose_user_name``
                      option set to ``True`` in the ``config/galaxy.ini``
                      configuration file. This parameter is silently ignored in
                      Galaxy ``release_15.10`` and earlier.
      --f_any TEXT    filter for user email or name. Each filter will be active for
                      non-admin users only if the Galaxy instance has the
                      corresponding ``expose_user_*`` option set to ``True`` in the
                      ``config/galaxy.ini`` configuration file. This parameter is
                      silently ignored in Galaxy ``release_15.10`` and earlier.
      -h, --help      Show this message and exit.
    

``show_user`` command
---------------------

**Usage**::

    parsec users show_user [OPTIONS] USER_ID

**Help**

Display information about a user.

**Options**::


      --deleted   whether to return results for a deleted user
      -h, --help  Show this message and exit.
    
