import click
from chakin.commands.analysis.add_analysis import cli as func0


@click.group()
def cli():
    """
    Access to the chado analysis table
    """
    pass


cli.add_command(func0)
