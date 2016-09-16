from os.path import dirname
import sys, pkgutil

for importer, package_name, _ in pkgutil.iter_modules([dirname(__file__)]):
    full_package_name = '%s.%s' % ('chado.actions', package_name)
    if full_package_name not in sys.modules:
        module = importer.find_module(package_name).load_module(full_package_name)
