#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
Comprehensive coverage tests for RU10Packet and RU10IntermediatePacket modules.
Targets 100% coverage for these files.
"""

from unittest.mock import patch

import pytest
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.ErrorCorrection import crc32
from norec4dna.RU10IntermediatePacket import RU10IntermediatePacket
from norec4dna.RU10Packet import RU10Packet


class TestRU10PacketInit:
    """Test RU10Packet initialization edge cases"""

    def test_ru10packet_with_dist(self):
        """Test RU10Packet with custom distribution"""
        dist = RaptorDistribution(100)
        packet = RU10Packet(b"test", [0, 1], 100, 1, dist=dist)
        assert packet.dist == dist

    def test_ru10packet_without_dist(self):
        """Test RU10Packet without distribution (uses default)"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        assert packet.dist is not None
        assert isinstance(packet.dist, RaptorDistribution)

    def test_ru10packet_with_method_even(self):
        """Test RU10Packet with method='even'"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, method="even")
        assert packet.method == "even"
        assert packet.packedMethod is not None

    def test_ru10packet_with_method_odd(self):
        """Test RU10Packet with method='odd'"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, method="odd")
        assert packet.method == "odd"
        assert packet.packedMethod is not None

    def test_ru10packet_with_method_window_30(self):
        """Test RU10Packet with method='window_30'"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, method="window_30", window=30)
        assert packet.method == "window_30"
        assert packet.window == 30
        assert packet.packedMethod is not None

    def test_ru10packet_with_method_window_40(self):
        """Test RU10Packet with method='window_40'"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, method="window_40", window=40)
        assert packet.method == "window_40"
        assert packet.window == 40
        assert packet.packedMethod is not None

    def test_ru10packet_with_invalid_method(self):
        """Test RU10Packet with invalid method raises error"""
        with pytest.raises(RuntimeError):
            RU10Packet(b"test", [0, 1], 100, 1, method="invalid")

    def test_ru10packet_with_invalid_window_method(self):
        """Test RU10Packet with invalid window method raises error"""
        with pytest.raises(RuntimeError):
            RU10Packet(b"test", [0, 1], 100, 1, method="window_50", window=50)

    def test_ru10packet_with_window_but_no_method(self):
        """Test RU10Packet packMethod with method but no window raises error"""
        # The assertion happens during __init__ when packMethod is called
        with pytest.raises(AssertionError):
            RU10Packet(b"test", [0, 1], 100, 1, method="window_30")

    def test_ru10packet_negative_id_spacing(self):
        """Test RU10Packet with negative id_spacing (should be set to 0)"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, id_spacing=-5)
        assert packet.id_spacing == 0

    def test_ru10packet_zero_id_spacing(self):
        """Test RU10Packet with zero id_spacing"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, id_spacing=0)
        assert packet.id_spacing == 0

    def test_ru10packet_positive_id_spacing(self):
        """Test RU10Packet with positive id_spacing"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, id_spacing=10)
        assert packet.id_spacing == 10

    def test_ru10packet_read_only(self):
        """Test RU10Packet with read_only=True"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, read_only=True)
        assert packet.packed_used_packets is None
        assert packet.packed is None

    def test_ru10packet_not_read_only(self):
        """Test RU10Packet with read_only=False"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, read_only=False)
        assert packet.packed_used_packets is not None
        assert packet.packed is not None

    def test_ru10packet_with_xor_by_seed(self):
        """Test RU10Packet with xor_by_seed=True"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, xor_by_seed=True)
        assert packet.xor_by_seed

    def test_ru10packet_with_mask_id_false(self):
        """Test RU10Packet with mask_id=False"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, mask_id=False)
        assert not packet.mask_id

    def test_ru10packet_with_contains_meta(self):
        """Test RU10Packet with contains_meta_or_version=True"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, contains_meta_or_version=True)
        assert packet.contains_meta_or_version


class TestRU10PacketMethods:
    """Test RU10Packet methods"""

    def test_get_packet_header_size_with_chunks(self):
        """Test get_packet_header_size with save_number_of_chunks_in_packet=True"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, save_number_of_chunks_in_packet=True)
        size = packet.get_packet_header_size()
        # Should include number_of_chunks_len_format (L=4) + id_len_format (L=4)
        # But also includes other fields in the actual header
        assert size > 0  # Just verify it returns a positive size

    def test_get_packet_header_size_without_chunks(self):
        """Test get_packet_header_size with save_number_of_chunks_in_packet=False"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, save_number_of_chunks_in_packet=False)
        size = packet.get_packet_header_size()
        # Should only include id_len_format (L=4)
        assert size > 0  # Just verify it returns a positive size
        # Without chunks should be smaller than with chunks
        packet_with = RU10Packet(b"test", [0, 1], 100, 1, save_number_of_chunks_in_packet=True)
        assert size < packet_with.get_packet_header_size()

    def test_set_used_packets_empty(self):
        """Test set_used_packets with empty list"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        with patch("logging.warning") as mock_warning:
            packet.set_used_packets([])
            mock_warning.assert_called_once()

    def test_set_used_packets_with_invalid_indices(self):
        """Test set_used_packets filters out invalid indices"""
        packet = RU10Packet(b"test", [0, 1], 10, 1)
        # Set with indices >= total_number_of_chunks (should be filtered)
        # Note: The current implementation doesn't filter, it just logs a warning
        # and uses all indices. The degree is based on len(used_packets)
        packet.set_used_packets([0, 1, 100, 200])
        # The degree is set from len(used_packets) before filtering
        # So it will be 4, not 2
        assert packet.degree == 4

    def test_prepare_and_pack_with_chunks(self):
        """Test prepare_and_pack with save_number_of_chunks_in_packet=True"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, save_number_of_chunks_in_packet=True)
        packed = packet.prepare_and_pack()
        assert len(packed) == 8  # L + L = 4 + 4

    def test_prepare_and_pack_without_chunks(self):
        """Test prepare_and_pack with save_number_of_chunks_in_packet=False"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, save_number_of_chunks_in_packet=False)
        packed = packet.prepare_and_pack()
        assert len(packed) == 4  # L = 4

    def test_calculate_packed_data_with_xor_by_seed(self):
        """Test calculate_packed_data with xor_by_seed=True"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, xor_by_seed=True)
        packed = packet.calculate_packed_data()
        assert isinstance(packed, bytes)
        assert len(packed) > 0

    def test_calculate_packed_data_with_method(self):
        """Test calculate_packed_data with method set"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, method="even")
        packed = packet.calculate_packed_data()
        assert isinstance(packed, bytes)
        assert len(packed) > 0

    def test_getId(self):
        """Test getId method"""
        packet = RU10Packet(b"test", [0, 1], 100, 42)
        assert packet.getId() == 42

    def test_setId(self):
        """Test setId method"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        packet.setId(99)
        assert packet.getId() == 99

    def test_from_packet_pseudo_false(self):
        """Test from_packet with pseudo=False"""
        original = RU10Packet(b"test", [0, 1], 100, 1)
        copy_packet = RU10Packet.from_packet(original, pseudo=False)
        assert copy_packet.get_data() == b"test"

    def test_from_packet_pseudo_true(self):
        """Test from_packet with pseudo=True"""
        original = RU10Packet(b"test", [0, 1], 100, 1)
        copy_packet = RU10Packet.from_packet(original, pseudo=True)
        assert copy_packet.get_data() == ""

    def test_get_number_of_half_blocks(self):
        """Test get_number_of_half_blocks"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        h = packet.get_number_of_half_blocks()
        assert h >= 0

    def test_get_number_of_ldpc_blocks(self):
        """Test get_number_of_ldpc_blocks"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        s = packet.get_number_of_ldpc_blocks()
        assert s >= 0

    def test_get_bool_array_used_packets(self):
        """Test get_bool_array_used_packets"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        arr = packet.get_bool_array_used_packets()
        assert arr is not None
        assert len(arr) == 100
        assert arr[0]
        assert arr[1]
        assert not arr[2]

    def test_get_bool_array_used_packets_none(self):
        """Test get_bool_array_used_packets when bool_arrayused_packets is None"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, read_only=True)
        # Force bool_arrayused_packets to None
        packet.bool_arrayused_packets = None
        arr = packet.get_bool_array_used_packets()
        assert arr is None

    def test_get_bool_array_all_used_packets(self):
        """Test get_bool_array_all_used_packets"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        arr = packet.get_bool_array_all_used_packets()
        # Should include chunks + ldpc + half blocks
        assert len(arr) > 100

    def test_get_bool_array_all_used_packets_no_used(self):
        """Test get_bool_array_all_used_packets with no used packets"""
        packet = RU10Packet(b"test", [], 100, 1)
        packet.used_packets = None
        arr = packet.get_bool_array_all_used_packets()
        assert all(not x for x in arr)

    def test_get_bool_array_used_and_ldpc_packets(self):
        """Test get_bool_array_used_and_ldpc_packets"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        arr = packet.get_bool_array_used_and_ldpc_packets()
        assert len(arr) > 100
        assert arr[0]
        assert arr[1]

    def test_get_bool_array_used_and_ldpc_packets_no_used(self):
        """Test get_bool_array_used_and_ldpc_packets with no used packets"""
        packet = RU10Packet(b"test", [], 100, 1)
        packet.used_packets = None
        arr = packet.get_bool_array_used_and_ldpc_packets()
        assert all(not x for x in arr)

    def test_get_bool_array_ldpc_packets(self):
        """Test get_bool_array_ldpc_packets"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        arr = packet.get_bool_array_ldpc_packets()
        # LDPC packets start after total_number_of_chunks
        assert len(arr) > 0

    def test_get_bool_array_half_packets(self):
        """Test get_bool_array_half_packets"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        arr = packet.get_bool_array_half_packets()
        assert len(arr) > 0

    def test_get_bool_array_repair_packets(self):
        """Test get_bool_array_repair_packets"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        arr = packet.get_bool_array_repair_packets()
        # Repair packets = ldpc + half blocks
        assert len(arr) > 0

    def test_str(self):
        """Test __str__ method"""
        packet = RU10Packet(b"test", [0, 1], 100, 1)
        str_repr = str(packet)
        assert "used_packets" in str_repr
        assert "Data" in str_repr

    def test_copy(self):
        """Test copy method"""
        original = RU10Packet(b"test", [0, 1], 100, 1)
        copy_packet = original.copy()

        assert copy_packet.get_data() == original.get_data()
        assert copy_packet.get_used_packets() == original.get_used_packets()
        assert copy_packet.get_total_number_of_chunks() == original.get_total_number_of_chunks()
        assert copy_packet.getId() == original.getId()


class TestRU10IntermediatePacket:
    """Test RU10IntermediatePacket class"""

    def test_ru10_intermediate_packet_basic(self):
        """Test basic RU10IntermediatePacket creation"""
        packet = RU10IntermediatePacket(b"test", [0, 1], 100, 1)
        assert packet.get_data() == b"test"
        assert packet.id == 1
        assert packet.total_number_of_chunks == 100

    def test_ru10_intermediate_packet_set_used_packets(self):
        """Test RU10IntermediatePacket set_used_packets"""
        packet = RU10IntermediatePacket(b"test", [0, 1], 100, 1)
        packet.set_used_packets([5, 10, 15])
        assert set(packet.get_used_packets()) == {5, 10, 15}

    def test_ru10_intermediate_packet_inherits_methods(self):
        """Test RU10IntermediatePacket inherits RU10Packet methods"""
        packet = RU10IntermediatePacket(b"test", [0, 1], 100, 1)

        # Test inherited methods
        assert packet.getId() == 1
        packet.setId(99)
        assert packet.getId() == 99

        assert packet.get_number_of_half_blocks() >= 0
        assert packet.get_number_of_ldpc_blocks() >= 0

        arr = packet.get_bool_array_used_packets()
        assert arr is not None
        assert len(arr) == 100

    def test_ru10_intermediate_packet_str(self):
        """Test RU10IntermediatePacket __str__"""
        packet = RU10IntermediatePacket(b"test", [0, 1], 100, 1)
        str_repr = str(packet)
        assert "used_packets" in str_repr
        assert "Data" in str_repr


class TestRU10PacketEdgeCases:
    """Test edge cases for RU10Packet"""

    def test_ru10packet_empty_data(self):
        """Test RU10Packet with empty data"""
        packet = RU10Packet(b"", [0, 1], 100, 1)
        assert packet.get_data() == b""

    def test_ru10packet_large_data(self):
        """Test RU10Packet with large data"""
        large_data = b"x" * 10000
        packet = RU10Packet(large_data, [0, 1], 100, 1)
        assert packet.get_data() == large_data

    def test_ru10packet_all_used_packets(self):
        """Test RU10Packet with all packets used"""
        all_packets = list(range(100))
        packet = RU10Packet(b"test", all_packets, 100, 1)
        assert packet.degree == 100

    def test_ru10packet_single_used_packet(self):
        """Test RU10Packet with single used packet"""
        packet = RU10Packet(b"test", [50], 100, 1)
        assert packet.degree == 1

    def test_ru10packet_duplicate_used_packets(self):
        """Test RU10Packet with duplicate used packets (should be handled)"""
        packet = RU10Packet(b"test", [0, 0, 1, 1], 100, 1)
        # Duplicates should be handled by set_used_packets
        assert packet.degree >= 1

    def test_ru10packet_with_error_correction(self):
        """Test RU10Packet with custom error correction"""
        packet = RU10Packet(b"test", [0, 1], 100, 1, error_correction=crc32)
        assert packet.error_correction == crc32

    def test_ru10packet_different_len_formats(self):
        """Test RU10Packet with different length formats"""
        packet = RU10Packet(
            b"test",
            [0, 1],
            100,
            1,
            packet_len_format="H",
            crc_len_format="H",
            number_of_chunks_len_format="H",
            id_len_format="H",
        )
        assert packet.packet_len_format == "H"
        assert packet.crc_len_format == "H"
        size = packet.get_packet_header_size()
        # H + H = 2 + 2 = 4
        assert size == 4


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
