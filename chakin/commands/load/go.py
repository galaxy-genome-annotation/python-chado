import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, dict_output


@click.command('go')
@click.argument("input", type=str)
@click.argument("organism_id", type=int)
@click.argument("analysis_id", type=int)
@click.option(
    "--query_type",
    help="The feature type (e.g. 'gene', 'mRNA', 'polypeptide', 'contig') of the query. It must be a valid Sequence Ontology term.",
    default="polypeptide",
    show_default=True,
    type=str
)
@click.option(
    "--match_on_name",
    help="Match features using their name instead of their uniquename",
    is_flag=True
)
@click.option(
    "--name_column",
    help="Column containing the feature identifiers (2, 3, 10 or 11; default=2).",
    default="2",
    show_default=True,
    type=int
)
@click.option(
    "--go_column",
    help="Column containing the GO id (default=5).",
    default="5",
    show_default=True,
    type=int
)
@click.option(
    "--re_name",
    help="Regular expression to extract the feature name from the input file (first capturing group will be used).",
    type=str
)
@click.option(
    "--skip_missing",
    help="Skip lines with unknown features or GO id instead of aborting everything.",
    is_flag=True
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, input, organism_id, analysis_id, query_type="polypeptide", match_on_name=False, name_column=2, go_column=5, re_name="", skip_missing=False):
    """Load GO annotation from a tabular file, in the same way as does the tripal_analysis_go module

Output:

    Number of inserted GO terms
    """
    return ctx.gi.load.go(input, organism_id, analysis_id, query_type=query_type, match_on_name=match_on_name, name_column=name_column, go_column=go_column, re_name=re_name, skip_missing=skip_missing)
