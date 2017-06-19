from __future__ import absolute_import
import os
import sys
import click
import json

from .io import error
from .config import read_global_config, global_config_path, set_global_config_path, get_instance  # noqa, ditto
from chakin import __version__  # noqa, ditto

CONTEXT_SETTINGS = dict(auto_envvar_prefix='CHAKIN', help_option_names=['-h', '--help'])


class Context(object):

    def __init__(self):
        self.verbose = False
        self.home = os.getcwd()
        self._global_config = None

    @property
    def global_config(self):
        if self._global_config is None:
            self._global_config = read_global_config()
        return self._global_config

    def log(self, msg, *args):
        """Logs a message to stderr."""
        if args:
            msg %= args
        click.echo(msg, file=sys.stderr)

    def vlog(self, msg, *args):
        """Logs a message to stderr only if verbose is enabled."""
        if self.verbose:
            self.log(msg, *args)

    def exit(self, exit_code):
        self.vlog("Exiting chakin with exit code [%d]" % exit_code)
        sys.exit(exit_code)


pass_context = click.make_pass_decorator(Context, ensure=True)
cmd_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          'commands'))


def list_cmds():
    rv = []
    for filename in os.listdir(cmd_folder):
        if filename.endswith('.py') and \
           filename.startswith('cmd_'):
            rv.append(filename[len("cmd_"):-len(".py")])
    rv.sort()
    return rv


def list_subcmds(parent):
    rv = []
    for filename in os.listdir(os.path.join(cmd_folder, parent)):
        if filename.endswith('.py') and \
           not filename.startswith('__'):
            rv.append(filename[:-len(".py")])
    rv.sort()
    return rv


def name_to_command(parent, name):
    try:
        if sys.version_info[0] == 2:
            if parent:
                parent = parent.encode('ascii', 'replace')
            name = name.encode('ascii', 'replace')

        if parent:
            mod_name = 'chakin.commands.%s.%s' % (parent, name)
        else:
            mod_name = 'chakin.commands.cmd_' + name
        mod = __import__(mod_name, None, None, ['cli'])
    except ImportError as e:
        error("Problem loading command %s, exception %s" % (name, e))
        return
    return mod.cli


class chakinCLI(click.MultiCommand):

    def list_commands(self, ctx):
        # We pre-calculate this so it works more nicely within packaged
        # versions of chakin. Please feel free to fix this?

        commands = ['init', 'organism', 'export', 'util', 'analysis']
        return commands

    def get_command(self, ctx, name):
        return name_to_command(None, name)


@click.command(cls=chakinCLI, context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__)
@click.option('-v', '--verbose', is_flag=True,
              help='Enables verbose mode.')
@click.option(
    "-i",
    "--instance",
    help='Name of instance in %s. This parameter can also be set via the environment variable CHADO_INSTANCE' % global_config_path(),
    default='__default',
    show_default=True,
    required=True
)
@click.option(
    "--path", "-f",
    help="config file path",
    type=str
)
@pass_context
def chakin(ctx, instance, verbose, path=None):
    """Command line wrappers around BioBlend functions. While this sounds
    unexciting, with chakin and jq you can easily build powerful command line
    scripts."""
    # set config_path if provided
    if path is not None and len(path) > 0:
        set_global_config_path(path)
    # We abuse this, knowing that calls to one will fail.
    try:
        ctx.gi = get_instance(instance)
    except TypeError:
        pass

    ctx.verbose = verbose


def json_loads(data):
    """Load json data, allowing - to represent stdin."""
    if data is None:
        return ""

    if data == "-":
        return json.load(sys.stdin)
    elif os.path.exists(data):
        with open(data, 'r') as handle:
            return json.load(handle)
    else:
        return json.loads(data)
