import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('export_gbk')
@click.argument("organism_id", type=int)
@pass_context
@custom_exception
@None_output
def cli(ctx, organism_id):
    """Export organism features as genbank

Output:

    None
    """
    return ctx.gi.export.export_gbk(organism_id)
