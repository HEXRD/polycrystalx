import sys
import os
from pathlib import Path
from importlib import import_module


def get_input_module(imod):
    """Process input string to get module

    Parameters
    ----------
    imod: str or Path
       name of input module

    RETURNS
    -------
    module
       input module specified by user
    """
    # Expect the user to run from directory containing the input module.
    sys.path.insert(0, str(Path.cwd()))

    # Allow user to specify module by path.
    if os.path.exists(imod):
        if imod.endswith("/"):
            imod = imod[:-1]
        path, ext = os.path.splitext(imod)
        imod = path.replace("/", ".")

    user_inputs = import_module(imod)

    return user_inputs
