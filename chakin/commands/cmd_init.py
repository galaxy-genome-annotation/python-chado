import os

import click

from chado import ChadoInstance
from chakin.cli import pass_context
from chakin import config
from chakin.io import warn, info

CONFIG_TEMPLATE = """## Chado's chakin: Global Configuration File.
# Each stanza should contian a single galaxy server to control.
#
# You can set the key __default to the name of a default instance
__default: local

local:
    dbhost: "%(dbhost)s"
    dbname: "%(dbname)s"
    dbuser: "%(dbuser)s"
    dbpass: "%(dbpass)s"
    dbport: "%(dbport)s"
    dbschema: "%(dbschema)s"
    debug: "%(debug)s"
"""

SUCCESS_MESSAGE = (
    "Ready to go! Type `chakin` to get a list of commands you can execute."
)


@click.command("config_init")
@pass_context
def cli(ctx, url=None, api_key=None, admin=False, **kwds):
    """Help initialize global configuration (in home directory)
    """
    # TODO: prompt for values someday.
    click.echo("""Welcome to Chado's Chakin!""")
    if os.path.exists(config.global_config_path()):
        info("Your chakin configuration already exists. Please edit it instead: %s" % config.global_config_path())
        return 0

    while True:
        apollo_url= click.prompt("Please entry your Apollo's URL")
        apollo_username = click.prompt("Please entry your Apollo Username")
        apollo_password = click.prompt("Please entry your Apollo Password", hide_input=True)
        info("Testing connection...")
        try:
            ai = ApolloInstance(apollo_url, apollo_username, apollo_password)
            try:
                ai.metrics.get_metrics()
                # Ok, success
                info("Ok! Everything looks good.")
                break
            except Exception as e:
                print(e)
                warn("Error, we could not access the configuration data for your instance.")
                should_break = click.prompt("Continue despite inability to contact this Apollo Instance? [y/n]")
                if should_break in ('Y', 'y'):
                    break
        except Exception as e:
            warn("Error, we could not access the configuration data for your instance.")
            should_break = click.prompt("Continue despite inability to contact this Apollo Instance? [y/n]")
            if should_break in ('Y', 'y'):
                break

    config_path = config.global_config_path()
    if os.path.exists(config_path):
        warn("File %s already exists, refusing to overwrite." % config_path)
        return -1
    with open(config_path, "w") as f:
        f.write(
            CONFIG_TEMPLATE % {
                'url': apollo_url,
                'username': apollo_username,
                'password': apollo_password,
            })
        info(SUCCESS_MESSAGE)
