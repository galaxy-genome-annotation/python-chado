import click
from chakin.commands.analysis.add_analysis import cli as add_analysis
from chakin.commands.analysis.delete_analyses import cli as delete_analyses
from chakin.commands.analysis.get_analyses import cli as get_analyses


@click.group()
def cli():
    """
    Access to the chado analysis table
    """
    pass


cli.add_command(add_analysis)
cli.add_command(delete_analyses)
cli.add_command(get_analyses)
