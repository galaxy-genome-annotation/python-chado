import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('gene_order')
@click.option(
    "--nuke",
    help="Removes all previous gene ordering data",
    is_flag=True
)
@pass_context
@custom_exception
@None_output
def cli(ctx, nuke=False):
    """Orders all the genes in the database by their order on their respective chromosomes in the gene_order table (for use in https://github.com/legumeinfo/lis_context_viewer/).

Output:

    None
    """
    return ctx.gi.phylogeny.gene_order(nuke=nuke)
