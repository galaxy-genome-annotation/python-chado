import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, list_output


@click.command('get_feature_cvterms')
@click.argument("feature_id", type=int)
@pass_context
@custom_exception
@list_output
def cli(ctx, feature_id):
    """Get cvterms associated with a feature

Output:

    Feature cvterms
    """
    return ctx.gi.feature.get_feature_cvterms(feature_id)
