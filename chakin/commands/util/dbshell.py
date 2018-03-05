import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('dbshell')
@pass_context
@custom_exception
@None_output
def cli(ctx):
    """Open a psql session to the database

Output:

    None
    """
    return ctx.gi.util.dbshell()
