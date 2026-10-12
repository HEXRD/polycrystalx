"""Loader of options inputs"""
import numpy as np

from ..inputs import options


class OptionsLoader:
    """Loader of options inputs

    Parameters
    ----------
    opts_inp: options.Options | None
        options input instance
    """

    def __init__(self, opts_inp):
        self.opts_inp = opts_inp
        # Get default options dictionary and update from loaded options.
        opts_d = options.Options()._asdict()
        if opts_inp is not None:
            new_d = opts_inp._asdict()
            opts_d.update(new_d)

        self.options = options.Options(**opts_d)
