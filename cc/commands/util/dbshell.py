import click
from cc.cli import pass_context, json_loads
from cc.decorators import chado_exception, None_output, _arg_split

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
