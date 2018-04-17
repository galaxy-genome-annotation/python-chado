from chakin.config import get_instance

ci = get_instance('local')


def setup_package():
    global ci
