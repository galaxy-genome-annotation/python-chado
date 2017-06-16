import click
from cc.commands.util.dbshell import cli as func0

@click.group()
def cli():
	pass

cli.add_command(func0)
