import numpy as np

def rounded(value, error, extra=1):
    """
    Return the value and error both rounded to the error's first digit
    """

    if error == 0:
        return (0,0)

    digit = int(np.ceil(-np.log10(error))) + extra
    assert np.round(error, digit) != 0
    return np.round(value, digit), np.round(error, digit)
