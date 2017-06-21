import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import chado_exception, dict_output, _arg_split

@click.command('add_analysis')
@click.argument("name", type=str)
@click.argument("program", type=str)
@click.argument("programversion", type=str)
@click.argument("algorithm")
@click.argument("sourcename", type=str)
@click.argument("sourceversion")
@click.argument("sourceuri", type=str)

@click.option(
    "--date_executed",
    help="analysis date_executed (yyyy-mm-dd)",
    type=str
)
@click.option(
    "--algorithn",
    help="analysis algorithn",
    type=str
)
@click.option(
    "--sourcevesion",
    help="analysis sourcevesion",
    type=str
)

@pass_context
@chado_exception
@dict_output
def cli(ctx, name, program, programversion, algorithm, sourcename, sourceversion, sourceuri, date_executed="", algorithn=None, sourcevesion=None):
    """Create an analysis

Output:

     Analysis information
        
    """
    kwargs = {}
    if algorithn and len(algorithn) > 0:
        kwargs['algorithn'] = algorithn
    if sourcevesion and len(sourcevesion) > 0:
        kwargs['sourcevesion'] = sourcevesion

    return ctx.gi.analysis.add_analysis(name, program, programversion, algorithm, sourcename, sourceversion, sourceuri, date_executed=date_executed, **kwargs)
