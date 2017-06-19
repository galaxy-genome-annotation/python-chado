import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import chado_exception, None_output, _arg_split

@click.command('export_fasta')
@click.argument("organism_id", type=int)

@click.option(
    "--file",
    help="If true, write to files in CWD",
    is_flag=True
)

@pass_context
@chado_exception
@None_output
def cli(ctx, organism_id, file=False):
    """Export reference sequences as fasta.

Output:

     None
        
    """
    return ctx.gi.export.export_fasta(organism_id, file=file)
