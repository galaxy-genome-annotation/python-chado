import click
from chakin.commands.analysis.add_analysis import cli as func0
from chakin.commands.analysis.delete_analyses import cli as func1
from chakin.commands.analysis.get_analyses import cli as func2


@click.group()
def cli():
    """
    Access to the chado analysis table
    """
    pass


cli.add_command(func0)
cli.add_command(func1)
cli.add_command(func2)
