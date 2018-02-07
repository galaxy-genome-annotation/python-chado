import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('load_fasta')
@click.argument("fasta", type=str)
@click.argument("analysis", type=int)
@click.argument("organism", type=int)
@click.argument("sequence_type", type=str)
@pass_context
@custom_exception
@None_output
def cli(ctx, fasta, analysis, organism, sequence_type):
    """Load features from a fasta file

Output:

    None
    """
    return ctx.gi.feature.load_fasta(fasta, analysis, organism, sequence_type)
