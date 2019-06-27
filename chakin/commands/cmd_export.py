import click
from chakin.commands.export.export_fasta import cli as export_fasta
from chakin.commands.export.export_gbk import cli as export_gbk
from chakin.commands.export.export_gff3 import cli as export_gff3


@click.group()
def cli():
    """
    Export data from the chado database
    """
    pass


cli.add_command(export_fasta)
cli.add_command(export_gbk)
cli.add_command(export_gff3)
