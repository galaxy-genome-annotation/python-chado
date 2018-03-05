import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, dict_output


@click.command('load_fasta')
@click.argument("fasta", type=str)
@click.argument("organism_id", type=int)
@click.option(
    "--sequence_type",
    help="Sequence type",
    default="contig",
    show_default=True,
    type=str
)
@click.option(
    "--analysis_id",
    help="Analysis ID",
    type=int
)
@click.option(
    "--re_name",
    help="Regular expression to extract the feature name from the fasta sequence id (first capturing group will be used).",
    type=str
)
@click.option(
    "--re_uniquename",
    help="Regular expression to extract the feature name from the fasta sequence id (first capturing group will be used).",
    type=str
)
@click.option(
    "--match_on_name",
    help="Match existing features using their name instead of their uniquename",
    is_flag=True
)
@click.option(
    "--update",
    help="Update existing feature with new sequence instead of throwing an error",
    is_flag=True
)
@click.option(
    "--db",
    help="External database to cross reference to.",
    type=int
)
@click.option(
    "--re_db_accession",
    help="Regular expression to extract an external database accession from the fasta sequence id (first capturing group will be used).",
    type=str
)
@click.option(
    "--rel_type",
    help="Relation type to parent feature ('part_of' or 'derives_from').",
    type=str
)
@click.option(
    "--re_parent",
    help="Regular expression to extract parent uniquename from the fasta sequence id (first capturing group will be used).",
    type=str
)
@click.option(
    "--parent_type",
    help="Sequence type of the parent feature",
    type=str
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, fasta, organism_id, sequence_type="contig", analysis_id="", re_name="", re_uniquename="", match_on_name=False, update=False, db="", re_db_accession="", rel_type="", re_parent="", parent_type=""):
    """Load features from a fasta file

Output:

    Number of inserted sequences
    """
    return ctx.gi.feature.load_fasta(fasta, organism_id, sequence_type=sequence_type, analysis_id=analysis_id, re_name=re_name, re_uniquename=re_uniquename, match_on_name=match_on_name, update=update, db=db, re_db_accession=re_db_accession, rel_type=rel_type, re_parent=re_parent, parent_type=parent_type)
