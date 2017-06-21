import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, None_output, _arg_split

@click.command('export_gff3')
@click.argument("organism_id", type=int)

@click.option(
    "--file",
    help=""
)

@pass_context
@custom_exception
@None_output
def cli(ctx, organism_id, file=False):
    """Export organism features as GFF3

Output:

     None
        
    """
    return ctx.gi.export.export_gff3(organism_id, file=file)
