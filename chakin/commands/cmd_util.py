import click
from chakin.commands.util.dbshell import cli as func0
from chakin.commands.util.launch_docker_image import cli as func1

@click.group()
def cli():
    pass

cli.add_command(func0)
cli.add_command(func1)
