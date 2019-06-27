import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, dict_output


@click.command('blast')
@click.argument("analysis_id", type=int)
@click.argument("input", type=str)
@click.option(
    "--blastdb",
    help="Name of the database blasted against (must be in the Chado db table)",
    type=str
)
@click.option(
    "--blastdb_id",
    help="ID of the database blasted against (must be in the Chado db table)",
    type=str
)
@click.option(
    "--blast_parameters",
    help="Blast parameters used to produce these results",
    type=str
)
@click.option(
    "--query_re",
    help="The regular expression that can uniquely identify the query name. This parameters is required if the feature name is not the first word in the blast query name.",
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
def cli(ctx, analysis_id, input, blastdb="", blastdb_id="", blast_parameters="", query_re="", query_type="polypeptide", query_uniquename=False):
    """Load a blast analysis, in the same way as does the tripal_analysis_blast module

Output:

    Number of processed hits
    """
    return ctx.gi.load.blast(analysis_id, input, blastdb=blastdb, blastdb_id=blastdb_id, blast_parameters=blast_parameters, query_re=query_re, query_type=query_type, query_uniquename=query_uniquename)
