workflows
=========

This section is auto-generated from the help text for the arrow command
``workflows``.


``cancel_invocation`` command
-----------------------------

**Usage**::

    parsec workflows cancel_invocation [OPTIONS] WORKFLOW_ID INVOCATION_ID

**Help**

Cancel the scheduling of a workflow.

**Options**::


      -h, --help  Show this message and exit.
    

``delete_workflow`` command
---------------------------

**Usage**::

    parsec workflows delete_workflow [OPTIONS] WORKFLOW_ID

**Help**

Delete a workflow identified by `workflow_id`.

**Options**::


      -h, --help  Show this message and exit.
    

``export_workflow_dict`` command
--------------------------------

**Usage**::

    parsec workflows export_workflow_dict [OPTIONS] WORKFLOW_ID

**Help**

Exports a workflow.

**Options**::


      -h, --help  Show this message and exit.
    

``export_workflow_json`` command
--------------------------------

**Usage**::

    parsec workflows export_workflow_json [OPTIONS] WORKFLOW_ID

**Help**

.. deprecated:: 0.9.0 Use :meth:`export_workflow_dict` instead.

**Options**::


      -h, --help  Show this message and exit.
    

``export_workflow_to_local_path`` command
-----------------------------------------

**Usage**::

    parsec workflows export_workflow_to_local_path [OPTIONS] WORKFLOW_ID

**Help**

Exports a workflow in JSON format to a given local path.

**Options**::


      --use_default_filename  If the use_default_name parameter is True, the
                              exported file will be saved as file_local_path/Galaxy-
                              Workflow-%s.ga, where %s is the workflow name. If
                              use_default_name is False, file_local_path is assumed
                              to contain the full file path including filename.
                              [default: True]
      -h, --help              Show this message and exit.
    

``get_invocations`` command
---------------------------

**Usage**::

    parsec workflows get_invocations [OPTIONS] WORKFLOW_ID

**Help**

Get a list containing all the workflow invocations corresponding to the specified workflow.

**Options**::


      -h, --help  Show this message and exit.
    

``get_workflow_inputs`` command
-------------------------------

**Usage**::

    parsec workflows get_workflow_inputs [OPTIONS] WORKFLOW_ID LABEL

**Help**

Get a list of workflow input IDs that match the given label. If no input matches the given label, an empty list is returned.

**Options**::


      -h, --help  Show this message and exit.
    

``get_workflows`` command
-------------------------

**Usage**::

    parsec workflows get_workflows [OPTIONS]

**Help**

Get all workflows or filter the specific one(s) via the provided ``name`` or ``workflow_id``. Provide only one argument, ``name`` or ``workflow_id``, but not both.

**Options**::


      --workflow_id TEXT  Encoded workflow ID (incompatible with ``name``)
      --name TEXT         Filter by name of workflow (incompatible with
                          ``workflow_id``). If multiple names match the given name,
                          all the workflows matching the argument will be returned.
      --published         if ``True``, return also published workflows
      -h, --help          Show this message and exit.
    

``import_shared_workflow`` command
----------------------------------

**Usage**::

    parsec workflows import_shared_workflow [OPTIONS] WORKFLOW_ID

**Help**

Imports a new workflow from the shared published workflows.

**Options**::


      -h, --help  Show this message and exit.
    

``import_workflow_dict`` command
--------------------------------

**Usage**::

    parsec workflows import_workflow_dict [OPTIONS] WORKFLOW_DICT

**Help**

Imports a new workflow given a dictionary representing a previously exported workflow.

**Options**::


      -h, --help  Show this message and exit.
    

``import_workflow_from_local_path`` command
-------------------------------------------

**Usage**::

    parsec workflows import_workflow_from_local_path [OPTIONS]

**Help**

Imports a new workflow given the path to a file containing a previously exported workflow.

**Options**::


      -h, --help  Show this message and exit.
    

``import_workflow_json`` command
--------------------------------

**Usage**::

    parsec workflows import_workflow_json [OPTIONS] WORKFLOW_JSON

**Help**

.. deprecated:: 0.9.0 Use :meth:`import_workflow_dict` instead.

**Options**::


      -h, --help  Show this message and exit.
    

``invoke_workflow`` command
---------------------------

**Usage**::

    parsec workflows invoke_workflow [OPTIONS] WORKFLOW_ID

**Help**

Invoke the workflow identified by ``workflow_id``. This will cause a workflow to be scheduled and return an object describing the workflow invocation.

**Options**::


      --inputs TEXT                   A mapping of workflow inputs to datasets and
                                      dataset collections. The datasets source can
                                      be a LibraryDatasetDatasetAssociation
                                      (``ldda``), LibraryDataset (``ld``),
                                      HistoryDatasetAssociation (``hda``), or
                                      HistoryDatasetCollectionAssociation
                                      (``hdca``).
      --params TEXT                   A mapping of non-datasets tool parameters (see
                                      below)
      --history_id TEXT               The encoded history ID where to store the
                                      workflow output. Alternatively,
                                      ``history_name`` may be specified to create a
                                      new history.
      --history_name TEXT             Create a new history with the given name to
                                      store the workflow output. If both
                                      ``history_id`` and ``history_name`` are
                                      provided, ``history_name`` is ignored. If
                                      neither is specified, a new 'Unnamed history'
                                      is created.
      --import_inputs_to_history      If ``True``, used workflow inputs will be
                                      imported into the history. If ``False``, only
                                      workflow outputs will be visible in the given
                                      history.
      --replacement_params TEXT       pattern-based replacements for post-job
                                      actions (see below)
      --allow_tool_state_corrections  If True, allow Galaxy to fill in missing tool
                                      state when running workflows. This may be
                                      useful for workflows using tools that have
                                      changed over time or for workflows built
                                      outside of Galaxy with only a subset of inputs
                                      defined.
      -h, --help                      Show this message and exit.
    

``run_invocation_step_action`` command
--------------------------------------

**Usage**::

    parsec workflows run_invocation_step_action [OPTIONS] WORKFLOW_ID

**Help**

nature of this action and what is expected will vary based on the the type of workflow step (the only currently valid action is True/False for pause steps).

**Options**::


      -h, --help  Show this message and exit.
    

``run_workflow`` command
------------------------

**Usage**::

    parsec workflows run_workflow [OPTIONS] WORKFLOW_ID

**Help**

Run the workflow identified by ``workflow_id``.

**Options**::


      --dataset_map TEXT          A mapping of workflow inputs to datasets. The
                                  datasets source can be a
                                  LibraryDatasetDatasetAssociation (``ldda``),
                                  LibraryDataset (``ld``), or
                                  HistoryDatasetAssociation (``hda``). The map must
                                  be in the following format: ``{'<input>': {'id':
                                  <encoded dataset ID>, 'src': '[ldda, ld, hda]'}}``
                                  (e.g. ``{'23': {'id': '29beef4fadeed09f', 'src':
                                  'ld'}}``)
      --params TEXT               A mapping of non-datasets tool parameters (see
                                  below)
      --history_id TEXT           The encoded history ID where to store the workflow
                                  output. Alternatively, ``history_name`` may be
                                  specified to create a new history.
      --history_name TEXT         Create a new history with the given name to store
                                  the workflow output. If both ``history_id`` and
                                  ``history_name`` are provided, ``history_name`` is
                                  ignored. If neither is specified, a new 'Unnamed
                                  history' is created.
      --import_inputs_to_history  If ``True``, used workflow inputs will be imported
                                  into the history. If ``False``, only workflow
                                  outputs will be visible in the given history.
      --replacement_params TEXT   pattern-based replacements for post-job actions
                                  (see below)
      -h, --help                  Show this message and exit.
    

``show_invocation`` command
---------------------------

**Usage**::

    parsec workflows show_invocation [OPTIONS] WORKFLOW_ID INVOCATION_ID

**Help**

Get a workflow invocation object representing the scheduling of a workflow. This object may be sparse at first (missing inputs and invocation steps) and will become more populated as the workflow is actually scheduled.

**Options**::


      -h, --help  Show this message and exit.
    

``show_invocation_step`` command
--------------------------------

**Usage**::

    parsec workflows show_invocation_step [OPTIONS] WORKFLOW_ID INVOCATION_ID

**Help**

See the details of a particular workflow invocation step.

**Options**::


      -h, --help  Show this message and exit.
    

``show_workflow`` command
-------------------------

**Usage**::

    parsec workflows show_workflow [OPTIONS] WORKFLOW_ID

**Help**

Display information needed to run a workflow.

**Options**::


      -h, --help  Show this message and exit.
    
