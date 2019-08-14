import click
from chakin.commands.util.dbshell import cli as dbshell
from chakin.commands.util.launch_docker_image import cli as launch_docker_image


@click.group()
def cli():
    """
    Some chado utilities
    """
    pass


cli.add_command(dbshell)
cli.add_command(launch_docker_image)
