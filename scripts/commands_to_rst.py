#!/usr/bin/env python

import os
import sys
from string import Template
import sqlalchemy

from click.testing import CliRunner

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, project_dir)

from chakin.cli import list_cmds, list_subcmds
from chakin import cli

chakin_cli = cli.chakin
runner = CliRunner()

COMMAND_TEMPLATE = Template('''
``${subcommand}`` command
${module_underline}

${command_help}
''')

COMMANDS_TEMPLATE = """========
Commands
========

Chakin is a set of wrappers for accessing Chado. Each utility is implemented as
a subcommand of ``chakin``. This section of the documentation
describes these commands.

.. toctree::
   :maxdepth: 0
"""

command_doc_dir = os.path.join("docs", "commands")
commands = COMMANDS_TEMPLATE

for command in list_cmds():
    if command == 'init':
        # Skip documenting init because it's special
        continue

    commands += "\n   commands/%s.rst" % command
    parent_doc_handle = open(os.path.join(command_doc_dir, command + ".rst"), "w")
    parent_doc_handle.write('%s\n' % command)
    parent_doc_handle.write('%s\n' % ('=' * len(command)))
    parent_doc_handle.write(Template("""
This section is auto-generated from the help text for the chakin command
``${command}``.

""").safe_substitute(command=command))


    for subcommand in list_subcmds(command):
        sqlalchemy.orm.clear_mappers()

        command_obj = cli.name_to_command(command, subcommand)

        function = command_obj.callback
        raw_rst = function.__doc__

        def clean_rst_line(line, remove=4):
            if line.startswith(" " * remove):
                return line[remove:]
            else:
                return line
        clean_rst = "\n".join(map(clean_rst_line, raw_rst.split("\n")))
        if 'Output:' in clean_rst:
            output_rst = clean_rst[clean_rst.index('Output:') + len('Output:'):].lstrip('\n')
            clean_rst = clean_rst[0:clean_rst.index('Output:')]
            output_rst = "\n".join([clean_rst_line(x, remove=5) for x in raw_rst.split("\n")])
        else:
            output_rst = ""

        result = runner.invoke(chakin_cli, [command, subcommand, "--help"])
        print(result)
        output = result.output
        lines = output.split("\n")
        new_lines = []
        help_lines = False
        option_lines = False
        output_lines = False

        for line in lines:
            if line.startswith("Usage: "):
                new_lines.append("**Usage**::\n\n    %s" % line[len("Usage: "):])
                new_lines.append("\n**Help**\n")
                new_lines.append(clean_rst)
                help_lines = True
                option_lines = False
                output_lines = False
            elif line.startswith("Options:"):
                help_lines = False
                option_lines = True
                output_lines = False
                new_lines.append("**Options**::\n\n")
            elif line.strip().startswith("Output:"):
                help_lines = False
                option_lines = False
                output_lines = True
                new_lines.append("**Output**\n\n")
                # I thought we already did this??
                new_lines.append(output_rst[output_rst.index('Output:') + len('Output:'):].lstrip('\n'))
            elif option_lines:
                new_lines.append("    %s" % line)
            elif output_lines:
                pass
        text = COMMAND_TEMPLATE.safe_substitute(
            command=command,
            subcommand=subcommand,
            command_help="\n".join(new_lines),
            module_underline = "-" * (len(subcommand) + len('```` command'))
        )
        parent_doc_handle.write(text)

open(os.path.join("docs", "commands.rst"), "w").write(commands)
