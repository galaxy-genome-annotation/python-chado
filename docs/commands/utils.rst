utils
=====

This section is auto-generated from the help text for the arrow command
``utils``.


``cmp`` command
---------------

**Usage**::

    parsec utils cmp [OPTIONS] METHOD CMP_WITH

**Help**

comparison tool. Exits if the value read from stdin does not pass the comparison test with the specified value.

method is the comparison method. One of these (lt, gt, eq, ne) which will
trigger a numerical comparison (cast to floats), or one of teq, tneq which
will trigger a string comparison.

cmp_with is the value to compare against. E.g. '100' or 'test'

e.g.::

    echo '5' | parsec utils cmp lt 10 # exit 0

    echo '5' | parsec utils cmp lt 1 # exit 1


**Options**::


      -h, --help  Show this message and exit.
    

``wait_on_invocation`` command
------------------------------

**Usage**::

    parsec utils wait_on_invocation [OPTIONS] WORKFLOW_ID INVOCATION_ID

**Help**

Given a workflow and invocation id, wait until that invocation is
complete (or one or more steps have errored)

This will exit with the following error codes:

- 0: done successfully
- 1: running (if --exit_early)
- 2: failure
- 3: unknown

**Options**::


      --exit_early         Exit immediately after checking status rather than
                           sleeping indefinitely
      --backoff_min FLOAT  Minimum time to sleep between checks, in seconds.
      --backoff_max FLOAT  Maximum time to sleep between checks, in seconds
      -h, --help           Show this message and exit.
    

``xunit_xargs`` command
-----------------------

**Usage**::

    parsec utils xunit_xargs [OPTIONS] _

**Help**

xargs look-alike that wraps output calls as XUnit XML

e.g.::

    parsec histories get_histories |             jq '.[].id' -r |             head -n 3 |             parsec utils xunit_xargs parsec histories get_status \| jq .percent_complete

will fetch the first three histories mentioned, and then pass them to xargs
to run ``parsec histories get_status [history_id] | jq .percent_complete``. This will
in turn produce XUnit XML that can be used in Jenkins or similar systems::

    <?xml version="1.0" ?>
    <testsuites errors="0" failures="0" tests="3" time="1.5944418907165527">
            <testsuite errors="0" failures="0" name="Parsec XX" skipped="0" tests="3" time="1.5944418907165527">
                    <testcase classname="parsec.histories.get_status.769f01a3981796db_|.jq..percent_complete" name="parsec.histories.get_status.769f01a3981796db_" time="0.604831">
                            <system-out>100</system-out>
                    </testcase>
                    <testcase classname="parsec.histories.get_status.83fbc32772cb5fcf_|.jq..percent_complete" name="parsec.histories.get_status.83fbc32772cb5fcf_" time="0.483556">
                            <system-out>100</system-out>
                    </testcase>
                    <testcase classname="parsec.histories.get_status.90c9282cb8718062_|.jq..percent_complete" name="parsec.histories.get_status.90c9282cb8718062_" time="0.506056">
                            <system-out>97.82608695652173</system-out>
                    </testcase>
            </testsuite>
    </testsuites>

**Options**::


      -h, --help  Show this message and exit.
    
