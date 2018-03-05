import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('export_gff3')
@click.argument("organism_id", type=int)
@pass_context
@custom_exception
@None_output
def cli(ctx, organism_id):
    """Export organism features as GFF3

Output:

    None
    """
    return ctx.gi.export.export_gff3(organism_id)
