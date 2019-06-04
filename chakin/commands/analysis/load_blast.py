import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, dict_output


@click.command('load_blast')
@click.argument("name", type=str)
@click.argument("program", type=str)
@click.argument("programversion", type=str)
@click.argument("sourcename", type=str)
@click.argument("blast_output", type=str)
@click.option(
    "--blast_ext",
    help="If looking for files in a directory, extension of the blast result files",
    type=str
)
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
    help="The feature type (e.g. 'gene', 'mRNA', 'contig') of the query. It must be a valid Sequence Ontology term.",
    type=str
)
@click.option(
    "--query_uniquename",
    help="Use this if the --query-re regular expression matches unique names instead of names in the database.",
    is_flag=True
)
@click.option(
    "--is_concat",
    help="If the blast result file is simply a list of concatenated blast results.",
    is_flag=True
)
@click.option(
    "--search_keywords",
    help="Extract keywords for Tripal search",
    is_flag=True
)
@click.option(
    "--no_parsed",
    help="Maximum number of hits to parse per feature. Default=all",
    default="all",
    show_default=True,
    type=str
)
@click.option(
    "--algorithm",
    help="analysis algorithm",
    type=str
)
@click.option(
    "--sourceversion",
    help="analysis sourceversion",
    type=str
)
@click.option(
    "--sourceuri",
    help="analysis sourceuri",
    type=str
)
@click.option(
    "--description",
    help="analysis description",
    type=str
)
@click.option(
    "--date_executed",
    help="analysis date_executed (yyyy-mm-dd)",
    type=str
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, name, program, programversion, sourcename, blast_output, blast_ext="", blastdb="", blastdb_id="", blast_parameters="", query_re="", query_type="", query_uniquename=False, is_concat=False, search_keywords=False, no_parsed="all", algorithm="", sourceversion="", sourceuri="", description="", date_executed=""):
    """Create an analysis

Output:

    Analysis information
    """
    return ctx.gi.analysis.load_blast(name, program, programversion, sourcename, blast_output, blast_ext=blast_ext, blastdb=blastdb, blastdb_id=blastdb_id, blast_parameters=blast_parameters, query_re=query_re, query_type=query_type, query_uniquename=query_uniquename, is_concat=is_concat, search_keywords=search_keywords, no_parsed=no_parsed, algorithm=algorithm, sourceversion=sourceversion, sourceuri=sourceuri, description=description, date_executed=date_executed)
