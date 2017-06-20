quotas
======

This section is auto-generated from the help text for the arrow command
``quotas``.


``create_quota`` command
------------------------

**Usage**::

    parsec quotas create_quota [OPTIONS] NAME DESCRIPTION AMOUNT OPERATION

**Help**

Create a new quota

**Options**::


      --default TEXT    Whether or not this is a default quota. Valid values are
                        ``no``, ``unregistered``, ``registered``. None is equivalent
                        to ``no``.  [default: no]
      --in_users TEXT   A list of user IDs or user emails.
      --in_groups TEXT  A list of group IDs or names.
      -h, --help        Show this message and exit.
    

``delete_quota`` command
------------------------

**Usage**::

    parsec quotas delete_quota [OPTIONS] QUOTA_ID

**Help**

Delete a quota

**Options**::


      -h, --help  Show this message and exit.
    

``get_quotas`` command
----------------------

**Usage**::

    parsec quotas get_quotas [OPTIONS]

**Help**

Get a list of quotas

**Options**::


      --deleted   Only return quota(s) that have been deleted
      -h, --help  Show this message and exit.
    

``show_quota`` command
----------------------

**Usage**::

    parsec quotas show_quota [OPTIONS] QUOTA_ID

**Help**

Display information on a quota

**Options**::


      --deleted   Search for quota in list of ones already marked as deleted
      -h, --help  Show this message and exit.
    

``undelete_quota`` command
--------------------------

**Usage**::

    parsec quotas undelete_quota [OPTIONS] QUOTA_ID

**Help**

Undelete a quota

**Options**::


      -h, --help  Show this message and exit.
    

``update_quota`` command
------------------------

**Usage**::

    parsec quotas update_quota [OPTIONS] QUOTA_ID

**Help**

Update an existing quota

**Options**::


      --name TEXT         Name for the new quota. This must be unique within a
                          Galaxy instance.
      --description TEXT  Quota description. If you supply this parameter, but not
                          the name, an error will be thrown.
      --amount TEXT       Quota size (E.g. ``10000MB``, ``99 gb``, ``0.2T``,
                          ``unlimited``)
      --operation TEXT    One of (``+``, ``-``, ``=``). If you wish to change this
                          value, you must also provide the ``amount``, otherwise it
                          will not take effect.
      --default TEXT      Whether or not this is a default quota. Valid values are
                          ``no``, ``unregistered``, ``registered``. Calling this
                          method with ``default="no"`` on a non-default quota will
                          throw an error. Not passing this parameter is equivalent
                          to passing ``no``.  [default: no]
      --in_users TEXT     A list of user IDs or user emails.
      --in_groups TEXT    A list of group IDs or names.
      -h, --help          Show this message and exit.
    
