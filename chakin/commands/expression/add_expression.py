import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, str_output


@click.command('add_expression')
@click.argument("organism_id", type=int)
@click.argument("analysis_id", type=int)
@click.argument("file_path", type=str)
@click.option(
    "--separator",
    help="Separating character in the matrix file (ex : ','). Default character is tab.",
    default="	",
    show_default=True,
    type=str
)
@click.option(
    "--unit",
    help="The units associated with the loaded values (ie, FPKM, RPKM, raw counts)",
    type=str
)
@click.option(
    "--query_type",
    help="The feature type (e.g. 'gene', 'mRNA', 'polypeptide', 'contig') of the query. It must be a valid Sequence Ontology term.",
    default="mRNA",
    show_default=True,
    type=str
)
@click.option(
    "--match_on_name",
    help="Match features using their name instead of their uniquename",
    is_flag=True
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
@str_output
def cli(ctx, organism_id, analysis_id, file_path, separator="	", unit="", query_type="mRNA", match_on_name=False, re_name="", skip_missing=False):
    """Add an expression matrix file to the database

Output:

    Number of expression data loaded
    """
    return ctx.gi.expression.add_expression(organism_id, analysis_id, file_path, separator=separator, unit=unit, query_type=query_type, match_on_name=match_on_name, re_name=re_name, skip_missing=skip_missing)
