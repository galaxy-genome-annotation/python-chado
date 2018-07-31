from chakin.config import get_instance

ci = get_instance('local')

ci_no_reflect = get_instance('local', no_reflect=True)


def setup_package():
    global ci
    global ci_no_reflect
