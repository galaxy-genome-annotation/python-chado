import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('load_gff')
@click.argument("gff", type=str)
@click.argument("analysis_id", type=int)
@click.argument("organism_id", type=int)
@click.option(
    "--landmark_type",
    help="Type of the landmarks (will speed up loading if provided, e.g. contig, should be a term of the Sequence ontology)",
    type=str
)
@click.option(
    "--re_protein",
    help="Replacement string for the protein name using capturing groups defined by --re_protein_capture",
    type=str
)
@click.option(
    "--re_protein_capture",
    help="Regular expression to capture groups in protein name to use in --re_protein (e.g. \"^(.*?)-R([A-Z]+)$\", default=\"^(.*?)$\")",
    default="^(.*?)$",
    show_default=True,
    type=str
)
@click.option(
    "--fasta",
    help="Path to a Fasta containing sequences for some features. When creating a feature, if its sequence is in this fasta file it will be loaded. Otherwise for mRNA and polypeptides it will be computed from the genome sequence (if available), otherwise it will be left empty.",
    type=str
)
@click.option(
    "--no_seq_compute",
    help="Disable the computation of mRNA and polypeptides sequences based on genome sequence and positions.",
    is_flag=True
)
@click.option(
    "--quiet",
    help="Hide progress information",
    is_flag=True
)
@click.option(
    "--add_only",
    help="Use this flag if you're not updating existing features, but just adding new features to the selected analysis and organism. It will speedup loading, and reduce memory usage, but might produce errors in case of already existing feature.",
    is_flag=True
)
@pass_context
@custom_exception
@None_output
def cli(ctx, gff, analysis_id, organism_id, landmark_type="", re_protein="", re_protein_capture="^(.*?)$", fasta="", no_seq_compute=False, quiet=False, add_only=False):
    """Load features from a gff file

Output:

    None
    """
    return ctx.gi.feature.load_gff(gff, analysis_id, organism_id, landmark_type=landmark_type, re_protein=re_protein, re_protein_capture=re_protein_capture, fasta=fasta, no_seq_compute=no_seq_compute, quiet=quiet, add_only=add_only)
