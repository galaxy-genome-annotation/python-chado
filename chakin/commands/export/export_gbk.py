import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import chado_exception, None_output, _arg_split

@click.command('export_gbk')
@click.argument("organism_id", type=int)


@pass_context
@chado_exception
@None_output
def cli(ctx, organism_id):
    """Export organism features as genbank

Output:

     None
        
    """
    return ctx.gi.export.export_gbk(organism_id)
