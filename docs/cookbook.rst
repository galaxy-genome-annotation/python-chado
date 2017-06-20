========
Cookbook
========

This page will contain more easy "recipes" for using  as time goes
on. Short tips and tricks that can help you use it more effectively, or
short recipes that can document how to do more complex tasks.

Talking to multiple Galaxies
----------------------------

If you are regularly switching between multiple Galaxy instances, you'll
probably want to take advantage of the environment variable for
specifying a Galaxy instance. E.g.:

.. code-block:: shell

    $ PARSEC_GALAXY_INSTANCE=uni-admin parsec config get_config | jq .brand
    "Internal"
    $ PARSEC_GALAXY_INSTANCE=uni-public parsec config get_config | jq .brand
    "Public"

You could easily set these at the top of a parsec script you've built
and all commands from there on would talk to the same Galaxy instance.

Capturing execution state as XUnit Output
-----------------------------------------

If you find yourself building a pipeline with ``parsec`` and ``jq``, you might find
yourself wanting to produce the output in a machine-legible format such as
XUnit. ``parsec`` now ships with a (very alpha) script to help with this.
``parsec utils xunit_xargs`` provides an ``xargs``-like experience, except it
produces XUnit formatted output. We'll run through a quick example of this:

.. code-block:: shell

    parsec histories get_histories | \
    jq '.[].id' -r | \
    head -n 3 | \
    parsec utils xunit_xargs parsec histories get_status \| jq .percent_complete \| parsec utils cmp eq 100

This command will fetch the first three histories, and then attempt to run
``parsec histories get_status <history_id> | jq .percent_complete`` for each
history id passed in.

.. code-block:: xml

    <?xml version="1.0" ?>
    <testsuites errors="0" failures="1" tests="3" time="1.6388022899627686">
      <testsuite errors="0" failures="1" name="Parsec XX" skipped="0" tests="3" time="1.6388022899627686">
        <testcase classname="parsec.histories.get_status.769f01a3981796db_|.jq..percent_complete.|.parsec.utils.cmp.eq.100" name="parsec.histories.get_status.769f01a3981796db_" time="0.537762"/>
        <testcase classname="parsec.histories.get_status.83fbc32772cb5fcf_|.jq..percent_complete.|.parsec.utils.cmp.eq.100" name="parsec.histories.get_status.83fbc32772cb5fcf_" time="0.534841"/>
        <testcase classname="parsec.histories.get_status.90c9282cb8718062_|.jq..percent_complete.|.parsec.utils.cmp.eq.100" name="parsec.histories.get_status.90c9282cb8718062_" time="0.566199">
          <failure message="Command 'parsec histories get_status 90c9282cb8718062 | jq .percent_complete | parsec utils cmp eq 100' returned non-zero exit status 1" type="failure">Traceback (most recent call last):

      File &quot;xunit_xargs.py&quot;, line 95, in cli
        output = check_output(' '.join(built_command), shell=True, stderr=stderr)

      File &quot;/usr/lib/python3.5/subprocess.py&quot;, line 626, in check_output
        **kwargs).stdout

      File &quot;/usr/lib/python3.5/subprocess.py&quot;, line 708, in run
        output=stdout, stderr=stderr)

    subprocess.CalledProcessError: Command 'parsec histories get_status 90c9282cb8718062 | jq .percent_complete | parsec utils cmp eq 100' returned non-zero exit status 1
    </failure>
          <system-err>97.82608695652173 != 100.0</system-err>
        </testcase>
      </testsuite>
    </testsuites>

Here we can see the example output, every history ID that went in came out as a
test case. One of them didn't pass a test we cared about and was marked as a
failure.
