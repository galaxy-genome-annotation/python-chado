import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, dict_output, _arg_split

@click.command('add_analysis')
@click.argument("name", type=str)
@click.argument("program", type=str)
@click.argument("programversion", type=str)
@click.argument("algorithm", type=str)
@click.argument("sourcename", type=str)
@click.argument("sourceversion", type=str)
@click.argument("sourceuri", type=str)

@click.option(
    "--date_executed",
    help="analysis date_executed (yyyy-mm-dd)",
    type=str
)

@pass_context
@custom_exception
@dict_output
def cli(ctx, name, program, programversion, algorithm, sourcename, sourceversion, sourceuri, date_executed=""):
    """Create an analysis

Output:

     Analysis information
        
    """
    return ctx.gi.analysis.add_analysis(name, program, programversion, algorithm, sourcename, sourceversion, sourceuri, date_executed=date_executed)
