requirements = [
    "lapy",
    "pyvista",
    "pyacvd",
    "nlopt==2.7.1"
]

# noinspection PyUnresolvedReferences
def ensure_requirements():
    # Imports are ensured at the top-level, once, when the module is initialized. So this function is a no-op,
    # but will explicitly fail if any requirement is missing. Note that this function does _not_ check
    # version constraints; use `install_requirements()` for that.

    import lapy
    import pyvista
    import pyacvd
    import nlopt


def install_requirements(upgrade=False):
    import slicer.util

    args = [*requirements]
    if upgrade:
        args.append("-U")

    slicer.util.pip_install(args)


try:
    ensure_requirements()
except ImportError:
    install_requirements()
