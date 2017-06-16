import click
from cc.commands.export.export_fasta import cli as func0
from cc.commands.export.export_gff3 import cli as func1

@click.group()
def cli():
	pass

cli.add_command(func0)
cli.add_command(func1)
