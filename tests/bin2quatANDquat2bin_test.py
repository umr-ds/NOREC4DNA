import pytest
from norec4dna.helper import bin2Quaternary, quaternary2Bin


@pytest.mark.parametrize("params", [("ACTG"), ("AAAA"), ("CCCC")])
def test_quat2bin_bin2quat(params):
    assert params == bin2Quaternary.byte2QUATS(
        int.from_bytes(quaternary2Bin.quats_to_bytes(params), "big")
    )
