import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import chado_exception, None_output, _arg_split

@click.command('dbshell')


@pass_context
@chado_exception
@None_output
def cli(ctx):
    """Open a psql session to the database

Output:

     None
        
    """
    return ctx.gi.util.dbshell()
