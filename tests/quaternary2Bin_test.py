import pytest
from contextlib import nullcontext as does_not_raise

from norec4dna.helper import quaternary2Bin


@pytest.mark.parametrize("params",
                         [('A', 0b00, does_not_raise()), ('C', 0b01, does_not_raise()), ('G', 0b10, does_not_raise()),
                          ('T', 0b11, does_not_raise()), ('X', None, pytest.raises(ValueError))])
def test_get_quarter_byte(params):
    with params[2]:
        assert quaternary2Bin.get_quarter_byte(params[0]) == params[1]
