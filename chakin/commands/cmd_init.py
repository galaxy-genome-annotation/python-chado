# coding: utf-8
import os

import click

from chado import ChadoInstance
from chakin.cli import pass_context
from chakin import config
from chakin.io import warn, info

CONFIG_TEMPLATE = """## Chado's chakin: Global Configuration File.
# Each stanza should contain a single chado server to control.
#
# You can set the key __default to the name of a default instance
__default: local

local:
    dbhost: "%(dbhost)s"
    dbname: "%(dbname)s"
    dbuser: "%(dbuser)s"
    dbpass: "%(dbpass)s"
    dbport: "%(dbport)s"
    dbschema: "%(schema)s"
"""

SUCCESS_MESSAGE = (
    "Ready to go! Type `chakin` to get a list of commands you can execute."
)


@click.command("config_init")
@pass_context
def cli(ctx, url=None, api_key=None, admin=False, **kwds):
    """Help initialize global configuration (in home directory)
    """

    click.echo("""Welcome to Chado's Chakin! (茶巾)""")
    if os.path.exists(config.global_config_path()):
        info("Your chakin configuration already exists. Please edit it instead: %s" % config.global_config_path())
        return 0

    while True:
        # Check environment
        dbhost = click.prompt("PGHOST")
        dbname = click.prompt("PGDATABASE")
        dbuser = click.prompt("PGUSER")
        dbpass = click.prompt("PGPASS", hide_input=True)
        dbport = click.prompt("PGPORT")
        schema = click.prompt("PGSCHEMA")

        info("Testing connection...")
        try:
            instance = ChadoInstance(dbhost=dbhost, dbname=dbname, dbuser=dbuser, dbpass=dbpass, dbport=dbport, dbschema=schema)
            # We do a connection test during startup.
            info("Ok! Everything looks good.")
            break
        except Exception as e:
            warn("Error, we could not access the configuration data for your instance: %s", e)
            should_break = click.prompt("Continue despite inability to contact this instance? [y/n]")
            if should_break in ('Y', 'y'):
                break

    config_path = config.global_config_path()
    if os.path.exists(config_path):
        warn("File %s already exists, refusing to overwrite." % config_path)
        return -1

    with open(config_path, "w") as f:
        f.write(CONFIG_TEMPLATE % {
            'dbhost': dbhost,
            'dbname': dbname,
            'dbuser': dbuser,
            'dbpass': dbpass,
            'dbport': dbport,
            'schema': schema,
        })
        info(SUCCESS_MESSAGE)
