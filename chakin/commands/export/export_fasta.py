import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('export_fasta')
@click.argument("organism_id", type=int)
@click.option(
    "--file",
    help="If true, write to files in CWD",
    is_flag=True
)
@pass_context
@custom_exception
@None_output
def cli(ctx, organism_id, file=False):
    """Export reference sequences as fasta.

Output:

    None
    """
    return ctx.gi.export.export_fasta(organism_id, file=file)
