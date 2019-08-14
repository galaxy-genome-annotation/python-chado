import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, dict_output


@click.command('interpro')
@click.argument("analysis_id", type=int)
@click.argument("input", type=str)
@click.option(
    "--parse_go",
    help="Load GO annotation to the database",
    is_flag=True
)
@click.option(
    "--query_re",
    help="The regular expression that can uniquely identify the query name. This parameter is required if the feature name is not the first word in the blast query name.",
    type=str
)
@click.option(
    "--query_type",
    help="The feature type (e.g. 'gene', 'mRNA', 'polypeptide', 'contig') of the query. It must be a valid Sequence Ontology term.",
    default="polypeptide",
    show_default=True,
    type=str
)
@click.option(
    "--query_uniquename",
    help="Use this if the --query-re regular expression matches unique names instead of names in the database.",
    is_flag=True
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, analysis_id, input, parse_go=False, query_re="", query_type="polypeptide", query_uniquename=False):
    """Load a blast analysis, in the same way as does the tripal_analysis_intepro module

Output:

    Number of processed hits
    """
    return ctx.gi.load.interpro(analysis_id, input, parse_go=parse_go, query_re=query_re, query_type=query_type, query_uniquename=query_uniquename)
