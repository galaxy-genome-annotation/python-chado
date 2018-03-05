import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, list_output


@click.command('add_cvterms')
@pass_context
@custom_exception
@list_output
def cli(ctx):
    """Make sure required cvterms are loaded

Output:

    created cvterms
    """
    return ctx.gi.phylogeny.add_cvterms()
