import click
from cc.commands.organism.delete_all_organisms import cli as func0
from cc.commands.organism.add_organism import cli as func1
from cc.commands.organism.get_organisms import cli as func2

@click.group()
def cli():
	pass

cli.add_command(func0)
cli.add_command(func1)
cli.add_command(func2)
