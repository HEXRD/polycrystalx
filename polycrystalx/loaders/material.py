"""Material input loader"""


class LinearElasticity:
    """Loader for LinearElasticity

    Parameters
    ----------
    userinput: module
       user input module
    """
    def __init__(self, userinput):
        self.materials = userinput.materials


# Slip has the Same interface as Elasticty.
Slip = LinearElasticity
