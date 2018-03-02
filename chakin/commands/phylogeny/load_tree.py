import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, dict_output


@click.command('load_tree')
@click.argument("newick", type=str)
@click.argument("analysis_id", type=int)
@click.option(
    "--name",
    help="The name given to the phylotree entry in the database (default=<filename>)",
    type=str
)
@click.option(
    "--xref_db",
    help="The name of the db to link dbxrefs for the trees (default: \"null\")",
    default="null",
    show_default=True,
    type=str
)
@click.option(
    "--xref_accession",
    help="The accession to use for dbxrefs for the trees (assumed same as name unless otherwise specified)",
    type=str
)
@click.option(
    "--match_on_name",
    help="Match polypeptide features using their name instead of their uniquename",
    is_flag=True
)
@click.option(
    "--prefix",
    help="Comma-separated list of prefix to be removed from identifiers (e.g species prefixes when using loading OrthoFinder output)",
    type=str
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, newick, analysis_id, name="", xref_db="null", xref_accession="", match_on_name=False, prefix=""):
    """Load a phylogenetic tree (Newick format) into Chado db

Output:

    Number of inserted trees
    """
    return ctx.gi.phylogeny.load_tree(newick, analysis_id, name=name, xref_db=xref_db, xref_accession=xref_accession, match_on_name=match_on_name, prefix=prefix)
