import click
from chakin.commands.organism.add_organism import cli as add_organism
from chakin.commands.organism.delete_all_organisms import cli as delete_all_organisms
from chakin.commands.organism.delete_organisms import cli as delete_organisms
from chakin.commands.organism.get_organisms import cli as get_organisms


@click.group()
def cli():
    """
    Access to the chado organism table
    """
    pass


cli.add_command(add_organism)
cli.add_command(delete_all_organisms)
cli.add_command(delete_organisms)
cli.add_command(get_organisms)
