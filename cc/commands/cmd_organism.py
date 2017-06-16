import click
from cc.commands.organism.get_organisms import cli as func0

@click.group()
def cli():
	pass

cli.add_command(func0)
