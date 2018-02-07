import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('load_gff')
@click.argument("gff", type=str)
@click.argument("analysis", type=int)
@click.argument("organism", type=int)
@click.argument("sequence_type", type=str)
@pass_context
@custom_exception
@None_output
def cli(ctx, gff, analysis, organism, sequence_type):
    """Load features from a gff file

Output:

    None
    """
    return ctx.gi.feature.load_gff(gff, analysis, organism, sequence_type)
