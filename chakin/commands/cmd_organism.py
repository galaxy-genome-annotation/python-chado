import click
from chakin.commands.organism.add_organism import cli as func0
from chakin.commands.organism.delete_all_organisms import cli as func1
from chakin.commands.organism.delete_organisms import cli as func2
from chakin.commands.organism.get_organisms import cli as func3


@click.group()
def cli():
    """
    Access to the chado organism table
    """
    pass


cli.add_command(func0)
cli.add_command(func1)
cli.add_command(func2)
cli.add_command(func3)
