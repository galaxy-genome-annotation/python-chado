import click
from chakin.commands.load.load_blast import cli as func0


@click.group()
def cli():
    pass


cli.add_command(func0)
