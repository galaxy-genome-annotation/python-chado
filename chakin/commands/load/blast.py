import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, dict_output


@click.command('blast')
@click.argument("analysis_id", type=int)
@click.argument("organism_id", type=int)
@click.argument("input", type=str)
@click.option(
    "--blastdb",
    help="Name of the database blasted against (must be in the Chado db table)",
    type=str
)
@click.option(
    "--blastdb_id",
    help="ID of the database blasted against (must be in the Chado db table)",
    type=int
)
@click.option(
    "--re_name",
    help="Regular expression to extract the feature name from the input file (first capturing group will be used).",
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
    "--match_on_name",
    help="Match features using their name instead of their uniquename",
    is_flag=True
)
@click.option(
    "--skip_missing",
    help="Skip lines with unknown features or GO id instead of aborting everything.",
    is_flag=True
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, analysis_id, organism_id, input, blastdb="", blastdb_id="", re_name="", query_type="polypeptide", match_on_name=False, skip_missing=False):
    """Load a blast analysis, in the same way as does the tripal_analysis_blast module

Output:

    Number of processed hits
    """
    return ctx.gi.load.blast(analysis_id, organism_id, input, blastdb=blastdb, blastdb_id=blastdb_id, re_name=re_name, query_type=query_type, match_on_name=match_on_name, skip_missing=skip_missing)
