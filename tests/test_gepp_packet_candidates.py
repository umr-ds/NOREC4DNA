import numpy as np
import pytest

from norec4dna.GEPP import GEPP


@pytest.mark.parametrize(
    "invalid, valid, multiple, expected",
    [
        ([0, 1], [], True, [True, True, True, False]),
        ([0], [], True, [True, True, False, False]),
        ([0, 1], [2], True, [True, False, True, False]),
        ([0, 1], [], False, [False, True, False, False]),
        ([0, 1], [2], False, [False, False, False, False]),
        ([], [], True, [False, False, False, False]),
        ([], [], False, [True, True, True, True]),
    ],
)
def test_packet_candidates(invalid, valid, multiple, expected):
    decoder = GEPP(np.eye(4, dtype=bool), np.zeros((4, 1), dtype=np.uint8))
    decoder.chunk_to_used_packets = np.array(
        [[1, 1, 0, 0], [0, 1, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]],
        dtype=bool,
    )
    original = decoder.chunk_to_used_packets.copy()

    actual = decoder.get_common_packets(invalid, valid, multi_error_packet_mode=multiple)

    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(decoder.chunk_to_used_packets, original)
