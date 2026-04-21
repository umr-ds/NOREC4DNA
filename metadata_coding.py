import argparse
import io
import logging
import random
import struct
import typing

import numpy
import numpy as np
from ConfigWorker import ConfigReadAndExecute
from norec4dna import RU10Decoder, RU10Encoder
from norec4dna.ErrorCorrection import get_error_correction_encode
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.helper.bin2Quaternary import string2QUATS
from norec4dna.helper.helper import calc_crc
from norec4dna.helper.quaternary2Bin import tranlate_quat_to_byte
from norec4dna.Packet import Packet
from norec4dna.RU10Packet import RU10Packet
from python_utils.types import deprecated

from .invivo_window_decoder import load_fasta

try:
    from norec4dna_multiversion.reconstruction import SemiAutomaticReconstructionToolkit
except ImportError:
    from semi_automatic_reconstruction_toolkit import (
        SemiAutomaticReconstructionToolkit,  # type: ignore[no-redef]
    )


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

PoolLoadResult = typing.Tuple[
    ConfigReadAndExecute,
    RU10Decoder,
    SemiAutomaticReconstructionToolkit,
]


def get_padding_area_offset(
    header_chunk: HeaderChunk,
    filename_is_free_area: bool = False,
) -> typing.Tuple[int, int]:
    """
    Returns the position and length of the unused padding area from the header chunk
    @header_chunk: the header chunk to extract from
    @filename_is_free_area:
    @returns: Tuple containing the start offset of the padding area and the padding_length
    """
    padding_offset = struct.calcsize(
        header_chunk.last_chunk_len_format + header_chunk.last_chunk_len_format
    )
    if filename_is_free_area:
        padding_offset += (
            len(header_chunk.get_file_name()) - 1
        )  # -1 as we must ensure 1 byte contains 0x00 (zero terminated string!)
    padding_length = len(header_chunk.data) - padding_offset
    return padding_offset, padding_length


# Header should contain:
#   - # changed chunks
#   - hash of new file
#   - reference (e.g. by seed) to an artificial packet containing all packets (by seed) that were changed -> these would be inserted FIRST during decoding! (optional) -> this concept could also be used for inserting DIFF information (increasing and decreasing the size)


def update_header(
    encoder: RU10Encoder,
    semiautomatic_solver: SemiAutomaticReconstructionToolkit,
    new_hash: typing.Optional[bytes],
    new_filename: typing.Optional[str],
    new_version: typing.Optional[int],
) -> None:
    # new_header = HeaderChunk(packet)
    if new_filename is None:
        new_filename = ""

    old_header = HeaderChunk(
        Packet(
            semiautomatic_solver.decoder.GEPP.b[
                semiautomatic_solver.decoder.GEPP.result_mapping[0]
            ],
            {0},
            semiautomatic_solver.decoder.number_of_chunks,
            read_only=True,
        ),
        last_chunk_len_format="I",
        checksum_len_format=semiautomatic_solver.decoder.checksum_len_str,
    )
    new_header: HeaderChunk = old_header  # TODO: create new header with updated values
    tmp = np.frombuffer(encoder.chunks[0][1:], dtype=np.uint8).reshape(-1)
    checksum = calc_crc(io.BytesIO(tmp[: len(tmp) - old_header.last_chunk_length]))  # TODO!
    additional_payload = (
        b""  # TODO: add metadata or special information and pad it to fill the chunk.
    )
    old_header.update_header(new_filename, checksum, additional_payload)
    encoder.chunks.insert(0, np.frombuffer(new_header.data, dtype=np.uint8))
    pass


def calculate_new_last_chunk_packet(
    packet: RU10Packet, last_chunk_length: int, wanted_metadata: bytes
) -> RU10Packet:
    assert len(wanted_metadata) <= len(packet.data) - last_chunk_length
    # wanted_metadata = np.zeros_like(np_header_chunk)
    # metadata_array = np.zeros_like(np_header_chunk)
    start_idx = len(packet.data) - len(wanted_metadata)
    # metadata_array[start_idx:] = np.frombuffer(wanted_metadata, dtype=np.bool)
    np_packet = numpy.frombuffer(packet.data, dtype=np.uint8).copy()
    np_packet[start_idx:] = np.frombuffer(wanted_metadata, dtype=np.uint8)

    out_packet = RU10Packet(
        np_packet,
        used_packets=typing.cast(typing.Collection[int], packet.used_packets or set()),
        total_number_of_chunks=packet.total_number_of_chunks,
        id=packet.id,
        dist=packet.dist,
        read_only=False,
        error_correction=packet.error_correction,
        packet_len_format=packet.packet_len_format,
        crc_len_format=packet.crc_len_format,
        number_of_chunks_len_format=packet.number_of_chunks_len_format,
        id_len_format=packet.id_len_format,
        save_number_of_chunks_in_packet=packet.save_number_of_chunks_in_packet,
        method=packet.method if hasattr(packet, "method") else None,
        window=packet.window if hasattr(packet, "window") else None,
        prepend=packet.prepend,
        append=packet.append,
        xor_by_seed=packet.xor_by_seed,
        mask_id=packet.mask_id,
        id_spacing=packet.id_spacing,
    )
    return out_packet


def calculate_new_header_packet(
    packet: RU10Packet, header_chunk: HeaderChunk, wanted_metadata: bytes
) -> RU10Packet:
    padding_offset, padding_length = get_padding_area_offset(
        header_chunk, filename_is_free_area=True
    )
    assert len(wanted_metadata) <= padding_length
    np_header_chunk = numpy.frombuffer(header_chunk.data, dtype=np.uint8)
    # wanted_metadata = np.zeros_like(np_header_chunk)
    # metadata_array = np.zeros_like(np_header_chunk)
    start_idx = len(np_header_chunk) - len(wanted_metadata)
    # metadata_array[start_idx:] = np.frombuffer(wanted_metadata, dtype=np.bool)
    np_packet = numpy.frombuffer(packet.data, dtype=np.uint8).copy()
    np_packet[start_idx:] = np.frombuffer(wanted_metadata, dtype=np.uint8)

    out_packet = RU10Packet(
        np_packet,
        used_packets=typing.cast(typing.Collection[int], packet.used_packets or set()),
        total_number_of_chunks=packet.total_number_of_chunks,
        id=packet.id,
        dist=packet.dist,
        read_only=False,
        error_correction=packet.error_correction,
        packet_len_format=packet.packet_len_format,
        crc_len_format=packet.crc_len_format,
        number_of_chunks_len_format=packet.number_of_chunks_len_format,
        id_len_format=packet.id_len_format,
        save_number_of_chunks_in_packet=packet.save_number_of_chunks_in_packet,
        method=packet.method if hasattr(packet, "method") else None,
        window=packet.window if hasattr(packet, "window") else None,
        prepend=packet.prepend,
        append=packet.append,
        xor_by_seed=packet.xor_by_seed,
        mask_id=packet.mask_id,
        id_spacing=packet.id_spacing,
    )
    return out_packet


def create_new_metadata_packets(
    semiautomatic_solver: SemiAutomaticReconstructionToolkit,
    metadata: list[str],
    unwanted_metadata: list[str],
    replace_packets_with_header: bool = True,
    replace_packets_with_last_chunk: bool = False,
) -> typing.List[RU10Packet]:
    """
    Add metadata as short unique sequences to the DNA sequence pool. These can then be used to identify if the DNA pool contains files matching the searched metadata.
    The metadata can be e.g. perceptual hashes of images, keywords, author names, etc. using a predefined hashing algorithm.
    @semiautomatic_solver: SemiAutomaticReconstructionToolkit
    @metadata: list of DNA sequences to add as metadata must be shorter than the remaining available space in a chunk
    (padding area + filename area in the header chunk - 1 (zero termination of pruned filename)).
    Adding the same metadata multiple times will result in multiple packets containing the metadata increasing the resilience against errors.

    (filename field may be pruned to gain more space - one byte for 0x00 is required)
    1. Check if the metadata fits into the available space:
        a) only use padding area (right-aligned)
        b) use padding area + pruned filename area (right-aligned)
        c) raise an error
    2. Check if the metadata is unique (not already present in the existing payload)
    3. Calculate the diff between wanted and current are as "to_patch" and convert this to binary
    4. apply: tmp_header = old_header XOR <padding>diff<checksum_padding>
    5. new_packet = old_packet XOR tmp_header + (recalculate checksum / RS)
    """
    if not replace_packets_with_header and not replace_packets_with_last_chunk:
        raise RuntimeError(
            "At least one of replace_packets_with_header or replace_packets_with_last_chunk must be set!"
        )
    new_packets: typing.List[RU10Packet] = []
    packets_using_header: typing.List[RU10Packet] = []
    packets_using_last_chunk: typing.List[RU10Packet] = []
    if len(semiautomatic_solver.decoder.packets) == 0:
        raise ValueError("Decoder does not contain any packets!")
    tmp_gepp = semiautomatic_solver.decoder.GEPP.clone()
    tmp_gepp.solve()
    header_row = int(tmp_gepp.result_mapping[0][0])
    header_chunk = HeaderChunk(
        Packet(
            tmp_gepp.b[header_row],
            {0},
            semiautomatic_solver.decoder.number_of_chunks,
            read_only=True,
        ),
        last_chunk_len_format=semiautomatic_solver.last_chunk_len_format,
        checksum_len_format=semiautomatic_solver.decoder.checksum_len_str,
    )
    added_metadata_str: typing.Set[str] = set()
    for packet in semiautomatic_solver.decoder.packets:
        if type(packet) is str:
            # Raw DNA sequence as checksum / RS failed - indicator for metadata influence.. -> UNLESS we update the Checksum as well
            logger.warning("Invalid packet...")
            dna_str = packet
        else:  # RU10Packet or any other type:
            dna_str = string2QUATS(packet.data)
        for metadata_str in metadata:
            if metadata_str in dna_str:
                logger.warn(f"Found wanted metadata string {metadata_str} in a packet")

        for unwanted_metadata_str in unwanted_metadata:
            if unwanted_metadata_str in dna_str:
                logger.error(
                    f"Found unwanted metadata string {unwanted_metadata_str} in a packet! "
                    f"This will lead to false-positives during querying! - consider manually removing it"
                )
        if type(packet) is str:
            continue  # cannot determine chunk coverage for raw DNA strings
        packet = typing.cast(RU10Packet, packet)
        used_raw_chunks = semiautomatic_solver.decoder.removeAndXorAuxPackets(packet)
        if replace_packets_with_header and used_raw_chunks[0]:
            packets_using_header.append(packet)
        if replace_packets_with_last_chunk and used_raw_chunks[-1]:
            packets_using_last_chunk.append(packet)
    total_packets_to_use = packets_using_header + packets_using_last_chunk
    for i, metadata_str in enumerate(metadata):
        packet = total_packets_to_use[i % len(total_packets_to_use)]
        # select x packets if they contain the header
        if i < len(packets_using_header):
            new_packet = calculate_new_header_packet(
                packet, header_chunk, tranlate_quat_to_byte(metadata_str)
            )
        else:
            new_packet = calculate_new_last_chunk_packet(
                packet, header_chunk.last_chunk_length, tranlate_quat_to_byte(metadata_str)
            )
        new_packets.append(new_packet)
        added_metadata_str.add(metadata_str)

    return new_packets


@deprecated("Use encoder_from_decoder instead!")
def encoder_from_sart(
    semiautomatic_solver: SemiAutomaticReconstructionToolkit,
    config_worker: ConfigReadAndExecute,
    rules=None,
) -> RU10Encoder:
    """
    Generate an encoder instance from an existing decoder
    IMPORTANT: Before using the encoder, prepend the new header to the encoder.chunks (...insert(0,header)) and call
               encoder.generate_intermediate_blocks()
    """
    # file, number_of_chunks, distribution: Distribution, insert_header=True, pseudo_decoder=None,
    #             chunk_size=0, rules=None, error_correction=nocode, packet_len_format="I", crc_len_format="L",
    #             number_of_chunks_len_format="L", id_len_format="L", save_number_of_chunks_in_packet=True,
    #             mode_1_bmp=False, prepend="", append="", drop_upper_bound=1.0, keep_all_packets=False,
    #             checksum_len_str=None, xor_by_seed=False, mask_id=True, id_spacing=0):
    sctn_config = config_worker.config[config_worker.config.sections()[0]]
    chunk_size = sctn_config.getint("chunk_size")
    number_of_chunks_len_format = sctn_config.get("number_of_chunks_len_format", "I")
    id_len_format = sctn_config.get("id_len_format", "I")
    save_number_of_chunks_in_packet = sctn_config.getboolean("savenumberofchunks", False)
    xor_by_seed = sctn_config.getboolean("xor_by_seed", True)
    mask_id = sctn_config.getboolean("mask_id", True)
    number_of_chunks = semiautomatic_solver.decoder.number_of_chunks
    id_spacing = sctn_config.getint("id_spacing", 0)
    crc_len_format = sctn_config.get("crc_len_format", "L")
    checksum_len_str = sctn_config.get("checksum_len_str", "I")
    last_chunk_len_str = sctn_config.get("last_chunk_len_str", "I")
    semiautomatic_solver.decoder.solve()
    semiautomatic_solver.decoder.populate_header_chunk(last_chunk_len_str=last_chunk_len_str)
    dist = semiautomatic_solver.decoder.distribution
    assert dist is not None
    ec_name = sctn_config.get("error_correction", "nocode")
    repair_symbols = sctn_config.getint("repair_symbols", 2)
    encode_error_correction = get_error_correction_encode(ec_name, repair_symbols)
    ru10_encoder = RU10Encoder(
        semiautomatic_solver.decoder.headerChunk.file_name.decode(),
        number_of_chunks,
        dist,
        semiautomatic_solver.decoder.use_headerchunk,
        None,
        0,
        rules,
        encode_error_correction,
        "I",
        crc_len_format,
        number_of_chunks_len_format,
        id_len_format,
        save_number_of_chunks_in_packet,
        False,
        "",
        "",
        1.0,
        True,
        "",
        xor_by_seed,
        mask_id,
        id_spacing,
        last_chunk_len_str,
        repair_symbols,
    )
    ru10_encoder.chunk_size = (
        chunk_size  # avoid overwriting the number of chunks by postponing the chunk_size setup
    )
    ru10_encoder.checksum_len_str = checksum_len_str
    ru10_encoder.checksum = semiautomatic_solver.decoder.headerChunk.checksum
    # remove the header and any decoded rows after the last chunk. Padding will be auto-applied as we use GEPP.b
    ru10_encoder.chunks = [
        x for x in semiautomatic_solver.decoder.GEPP.b[0:number_of_chunks]
    ]  # TODO: check for off-by-one due to header
    ru10_encoder.generate_intermediate_blocks()  # missing command from the prepare block
    # (we MUST not call prepare() as it would add an additional header!)
    for packet in semiautomatic_solver.decoder.packets:
        packet.packed_used_packets = packet.prepare_and_pack()
        packet.packed = packet.calculate_packed_data()
    ru10_encoder.encodedPackets = set(semiautomatic_solver.decoder.packets.copy())

    # ru10_encoder.chunks.insert(0,header)
    # ru10_encoder.generate_intermediate_blocks()

    return ru10_encoder


def encoder_from_decoder(
    decoder: RU10Decoder, config_worker: ConfigReadAndExecute, rules=None
) -> RU10Encoder:
    """
    Generate an encoder instance from an existing decoder
    IMPORTANT: Before using the encoder, prepend the new header to the encoder.chunks (...insert(0,header)) and call
               encoder.generate_intermediate_blocks()
    """
    # file, number_of_chunks, distribution: Distribution, insert_header=True, pseudo_decoder=None,
    #             chunk_size=0, rules=None, error_correction=nocode, packet_len_format="I", crc_len_format="L",
    #             number_of_chunks_len_format="L", id_len_format="L", save_number_of_chunks_in_packet=True,
    #             mode_1_bmp=False, prepend="", append="", drop_upper_bound=1.0, keep_all_packets=False,
    #             checksum_len_str=None, xor_by_seed=False, mask_id=True, id_spacing=0):
    sctn_config = config_worker.config[config_worker.config.sections()[0]]
    chunk_size = sctn_config.getint("chunk_size")
    number_of_chunks_len_format = sctn_config.get("number_of_chunks_len_format", "I")
    id_len_format = sctn_config.get("id_len_format", "I")
    save_number_of_chunks_in_packet = sctn_config.getboolean("savenumberofchunks", False)
    xor_by_seed = sctn_config.getboolean("xor_by_seed", True)
    mask_id = sctn_config.getboolean("mask_id", True)
    number_of_chunks = decoder.number_of_chunks
    id_spacing = sctn_config.getint("id_spacing", 0)
    crc_len_format = sctn_config.get("crc_len_format", "L")
    checksum_len_str = sctn_config.get("checksum_len_str", "I")
    last_chunk_len_str = sctn_config.get("last_chunk_len_str", "I")
    decoder.solve()
    decoder.populate_header_chunk(last_chunk_len_str=last_chunk_len_str)
    dist = decoder.distribution
    assert dist is not None
    ec_name = sctn_config.get("error_correction", "nocode")
    repair_symbols = sctn_config.getint("repair_symbols", 2)
    encode_error_correction = get_error_correction_encode(ec_name, repair_symbols)
    ru10_encoder = RU10Encoder(
        decoder.headerChunk.file_name.decode(),
        number_of_chunks,
        dist,
        decoder.use_headerchunk,
        None,
        0,
        rules,
        encode_error_correction,
        "I",
        crc_len_format,
        number_of_chunks_len_format,
        id_len_format,
        save_number_of_chunks_in_packet,
        False,
        "",
        "",
        1.0,
        True,
        "",
        xor_by_seed,
        mask_id,
        id_spacing,
        last_chunk_len_str,
        repair_symbols,
    )
    ru10_encoder.chunk_size = (
        chunk_size  # avoid overwriting the number of chunks by postponing the chunk_size setup
    )
    ru10_encoder.checksum_len_str = checksum_len_str
    ru10_encoder.checksum = decoder.headerChunk.checksum
    # remove the header and any decoded rows after the last chunk. Padding will be auto-applied as we use GEPP.b
    ru10_encoder.chunks = [
        x for x in decoder.GEPP.b[0:number_of_chunks]
    ]  # TODO: check for off-by-one due to header
    ru10_encoder.generate_intermediate_blocks()  # missing command from the prepare block
    # (we MUST not call prepare() as it would add an additional header!)
    for packet in decoder.packets:
        packet.packed_used_packets = packet.prepare_and_pack()
        packet.packed = packet.calculate_packed_data()
    ru10_encoder.encodedPackets = set(decoder.packets.copy())

    # ru10_encoder.chunks.insert(0,header)
    # ru10_encoder.generate_intermediate_blocks()

    return ru10_encoder


"""
def is_inconsistent_les(semiautomatic_solver):
    rank_a = semiautomatic_solver.calculate_rank_A()
    rank_aug = semiautomatic_solver.calculate_rank_augmented_matrix()
    return rank_a != rank_aug
"""

"""
def main(semiautomatic_solver):
    # Encode:
    # Load the ini file
    # Decode all packets and find packets that use chunk 0 (header chunk)
    # For the last Calculate: Wanted DNA Sequence
    # Decode:
    # Load the ini file
    # decode using semi_automatic_reconstruction_toolkit
    # extract and parse all version (using inconsistent
    # save all versions + metadata into the filesystem
    pass

    header = semiautomatic_solver.decoder.headerChunk
    ru10_encoder = encoder_from_decoder(semiautomatic_solver, cfg_worker)
    update_header(ru10_encoder, semiautomatic_solver, "", "", "")
    # generate_packets(semiautomatic_solver, ru10_encoder, 10, None)
"""


def parse_metadata_file(file: str) -> typing.List[str]:
    output: typing.List[str] = []
    fasta = typing.cast(typing.Dict[str, str], load_fasta(file))
    for name, meta in fasta.items():
        name, amount = name.split("_")
        output += int(amount) * [meta]
    return output


def replace_packets(
    encoder: RU10Encoder,
    new_packet_list: typing.List[RU10Packet],
    keep_percentage: float = 0.0,
    keep_one_header: bool = False,
    keep_one_last_chunk: bool = False,
    decoder: typing.Optional[RU10Decoder] = None,
) -> None:
    """
    Replace matching packets (same seed) in *encoder* with the metadata-carrying *new_packet_list*.

    @encoder:            Encoder whose ``encodedPackets`` set is updated in-place.
    @new_packet_list:    New (modified) packets to insert.
    @keep_percentage:    Fraction (0.0–1.0) of the old matching packets to keep alongside the new
                         ones, increasing redundancy.  Default 0.0 removes all old matching packets.
    @keep_one_header:    When True, always keep at least one original packet that covers the header
                         chunk (chunk 0), so the filename is never lost entirely.  Requires
                         *decoder* to be provided.
    @keep_one_last_chunk: When True, always keep at least one original packet that covers the last
                          chunk.  Requires *decoder* to be provided.
    @decoder:            ``RU10Decoder`` instance used for ``removeAndXorAuxPackets`` lookups.
                         Required when *keep_one_header* or *keep_one_last_chunk* is True.
    """
    if (keep_one_header or keep_one_last_chunk) and decoder is None:
        raise ValueError(
            "'decoder' must be provided when keep_one_header or keep_one_last_chunk is True"
        )

    # TODO: we _MUST_ make sure that SOME of the original packets with the header are kept - otherwise the filename gets lost!
    # TODO: we _MUST_ make sure that each selected packet with the header can be used to decode WITHOUT any other modified packets - otherwise we may end up with a header that is constructed using two modified packets making it impossible or very hard to resolve the deltas
    # - or: calculate all deltas as we did for the automatic reconstruction to detect the unique error deltas and the possible packets to find the delta of each packet!
    num_packets_to_keep: int = int(np.ceil(len(new_packet_list) * keep_percentage))
    new_packets_seeds = [x.id for x in new_packet_list]
    packets_to_remove: typing.List[RU10Packet] = []
    # find all packets with same seed as the new packets
    for packet in encoder.encodedPackets:
        if packet.id in new_packets_seeds:
            packets_to_remove.append(packet)
    random.shuffle(packets_to_remove)
    num_packets_to_keep_offset = 0
    if keep_one_header:
        special_packet_to_keep = -1
        for packet_to_keep in range(num_packets_to_keep):
            if decoder.removeAndXorAuxPackets(packets_to_remove[packet_to_keep])[0]:
                special_packet_to_keep = packet_to_keep
                break
        for potential_packet_to_keep in range(num_packets_to_keep, len(packets_to_remove)):
            if decoder.removeAndXorAuxPackets(packets_to_remove[potential_packet_to_keep])[0]:
                special_packet_to_keep = potential_packet_to_keep
                break
        if special_packet_to_keep == -1:
            raise RuntimeError("No packet with header chunk found to keep!")
        # swap first with current packet_to_keep:
        packets_to_remove[0], packets_to_remove[special_packet_to_keep] = (
            packets_to_remove[special_packet_to_keep],
            packets_to_remove[0],
        )
        num_packets_to_keep_offset += 1

    if keep_one_last_chunk:
        special_packet_to_keep = -1
        for packet_to_keep in range(num_packets_to_keep):
            if decoder.removeAndXorAuxPackets(packets_to_remove[packet_to_keep])[-1]:
                special_packet_to_keep = packet_to_keep
                break
        for potential_packet_to_keep in range(num_packets_to_keep, len(packets_to_remove)):
            if decoder.removeAndXorAuxPackets(packets_to_remove[potential_packet_to_keep])[-1]:
                special_packet_to_keep = potential_packet_to_keep
                break
        if special_packet_to_keep == -1:
            raise RuntimeError("No packet with last chunk found to keep!")
        # swap second with current packet_to_keep:
        packets_to_remove[1], packets_to_remove[special_packet_to_keep] = (
            packets_to_remove[special_packet_to_keep],
            packets_to_remove[1],
        )
        num_packets_to_keep_offset += 1
    if num_packets_to_keep < num_packets_to_keep_offset:
        num_packets_to_keep = num_packets_to_keep_offset
    to_remove = set(packets_to_remove[int(num_packets_to_keep) :])
    to_remove_ids = {id(p) for p in to_remove}
    # Rebuild without relying on equality lookups during remove
    encoder.encodedPackets = {p for p in encoder.encodedPackets if id(p) not in to_remove_ids}
    skipped_packets = 0
    # add all new packets:
    for packet in new_packet_list:
        if packet in encoder.encodedPackets:
            logger.error(
                f"Packet already in encoded Packets - this packet will be skipped: {packet}"
            )
            skipped_packets += 1
        encoder.encodedPackets.add(packet)
    logger.warning(
        f"There were a total of {skipped_packets} packets containing duplicates. "
        f"The the actual amount of copies may be different from what you requested. "
        f"This usually indicates that there the requested number of copies is larger that the total number of packets containing the header chunk."
    )


def init_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="python -m NOREC4DNA.metadata_coding",
        description=(
            "Embed searchable DNA metadata into a NOREC4DNA pool, query a pool for metadata "
            "sequences, or inspect padding capacity."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Embed one metadata sequence (inline):
  python -m NOREC4DNA.metadata_coding embed \\
      --ini pool.ini --metadata "ACGTACGTACGTACGT" --output pool_meta.fasta

  # Embed multiple sequences from a FASTA file:
  python -m NOREC4DNA.metadata_coding embed \\
      --ini pool.ini --metadata_file meta.fasta --output pool_meta.fasta

  # Keep 20%% of original matching packets (higher redundancy):
  python -m NOREC4DNA.metadata_coding embed \\
      --ini pool.ini --metadata "ACGTACGT" --keep_percentage 0.2 --output out.fasta

  # Search a pool for a metadata sequence:
  python -m NOREC4DNA.metadata_coding query \\
      --fasta pool.fasta --metadata "ACGTACGT"

  # Show available padding capacity for a pool:
  python -m NOREC4DNA.metadata_coding info --ini pool.ini
""",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ── embed ────────────────────────────────────────────────────────────────
    embed_p = subparsers.add_parser(
        "embed",
        help="Embed DNA metadata strings into a FASTA pool.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    embed_p.add_argument(
        "--ini", required=True, metavar="INI", help="INI config file of the pool to modify."
    )
    embed_p.add_argument(
        "--output",
        metavar="FASTA",
        default=None,
        help="Output FASTA path. Omit to overwrite the pool in-place.",
    )

    meta_g = embed_p.add_mutually_exclusive_group(required=True)
    meta_g.add_argument(
        "--metadata",
        metavar="DNA[,DNA,...]",
        help="Comma-separated DNA metadata sequences to embed.",
    )
    meta_g.add_argument(
        "--metadata_file",
        metavar="FILE",
        help="FASTA file of metadata sequences.  Header format: >name_N where N "
        "is the number of copies to embed (e.g. >hash_3 embeds the sequence 3 times).",
    )

    unwanted_g = embed_p.add_mutually_exclusive_group()
    unwanted_g.add_argument(
        "--unwanted",
        metavar="DNA[,DNA,...]",
        help="Comma-separated sequences that must NOT appear in any packet "
        "(logged as errors if found).",
    )
    unwanted_g.add_argument(
        "--unwanted_file",
        metavar="FILE",
        help="FASTA file of unwanted sequences (same format as --metadata_file).",
    )

    embed_p.add_argument(
        "--keep_percentage",
        type=float,
        default=0.0,
        metavar="FRAC",
        help="Fraction [0.0–1.0] of original matching packets to keep alongside "
        "the new metadata packets for extra redundancy.  Default removes all "
        "original matching packets.",
    )
    embed_p.add_argument(
        "--keep_one_header",
        action="store_true",
        help="Always preserve at least one original packet covering the header "
        "chunk (chunk 0) so the filename is never lost.",
    )
    embed_p.add_argument(
        "--keep_one_last_chunk",
        action="store_true",
        help="Always preserve at least one original packet covering the last chunk.",
    )
    embed_p.add_argument(
        "--use_last_chunk",
        action="store_true",
        help="Use last-chunk packets as carriers in addition to header-chunk packets.",
    )

    # ── query ────────────────────────────────────────────────────────────────
    query_p = subparsers.add_parser(
        "query",
        help="Search a FASTA pool for one or more metadata sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    query_p.add_argument("--fasta", required=True, metavar="FASTA", help="FASTA file to search.")

    qmeta_g = query_p.add_mutually_exclusive_group(required=True)
    qmeta_g.add_argument(
        "--metadata",
        metavar="DNA[,DNA,...]",
        help="Comma-separated metadata sequences to search for.",
    )
    qmeta_g.add_argument(
        "--metadata_file", metavar="FILE", help="FASTA file of metadata sequences to search for."
    )

    query_p.add_argument(
        "--count_only",
        action="store_true",
        help="Print only the match counts, not the matching sequence IDs.",
    )

    # ── info ─────────────────────────────────────────────────────────────────
    info_p = subparsers.add_parser(
        "info",
        help="Show padding capacity and pool statistics for metadata embedding.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    info_p.add_argument(
        "--ini", required=True, metavar="INI", help="INI config file of the pool to inspect."
    )

    return parser.parse_args()


# ── subcommand implementations ─────────────────────────────────────────────


def _load_pool(ini_file: str) -> PoolLoadResult:
    """Load INI, return (cfg_worker, decoder, semiautomatic_solver)."""
    cfg_worker = ConfigReadAndExecute(ini_file)
    decoder = typing.cast(RU10Decoder, cfg_worker.execute(return_decoder=True, skip_solve=True)[0])
    semiautomatic_solver = SemiAutomaticReconstructionToolkit(decoder)
    return cfg_worker, decoder, semiautomatic_solver


def _resolve_metadata(
    args_metadata: typing.Optional[str],
    args_metadata_file: typing.Optional[str],
) -> typing.List[str]:
    if args_metadata_file is not None:
        return parse_metadata_file(args_metadata_file)
    if args_metadata:
        return [s.strip() for s in args_metadata.split(",") if s.strip()]
    return []


def _cmd_embed(args: argparse.Namespace) -> int:
    """Embed metadata into a pool and save the result."""
    cfg_worker, decoder, semiautomatic_solver = _load_pool(args.ini)

    wanted = _resolve_metadata(
        getattr(args, "metadata", None), getattr(args, "metadata_file", None)
    )
    if not wanted:
        print("ERROR: no metadata sequences provided.", flush=True)
        return 1

    unwanted = _resolve_metadata(
        getattr(args, "unwanted", None), getattr(args, "unwanted_file", None)
    )

    use_header = True
    use_last = getattr(args, "use_last_chunk", False)
    if not use_last:
        # keep_one_last_chunk only makes sense when we also embed into last-chunk packets
        use_last = getattr(args, "keep_one_last_chunk", False)

    print(f"Embedding {len(wanted)} metadata sequence(s) into '{args.ini}' …", flush=True)

    encoder = encoder_from_decoder(decoder, cfg_worker)

    new_packets = create_new_metadata_packets(
        semiautomatic_solver,
        wanted,
        unwanted,
        replace_packets_with_header=use_header,
        replace_packets_with_last_chunk=use_last,
    )

    # Re-pack each modified packet so .packed reflects the new .data
    for pkt in new_packets:
        pkt.packed_used_packets = pkt.prepare_and_pack()
        pkt.packed = pkt.calculate_packed_data()

    replace_packets(
        encoder,
        new_packets,
        keep_percentage=args.keep_percentage,
        keep_one_header=args.keep_one_header,
        keep_one_last_chunk=args.keep_one_last_chunk,
        decoder=decoder if (args.keep_one_header or args.keep_one_last_chunk) else None,
    )

    # Determine output path
    import configparser

    ini_path = args.ini
    if args.output:
        out_fasta = args.output
    else:
        cfg = configparser.ConfigParser()
        cfg.read(ini_path)
        section = cfg.sections()[0]
        out_fasta = cfg[section].get("filename", section)

    encoder.save_packets_fasta(out_fasta, "", False)

    # Write companion INI — save_config_file always writes to {encoder.file}_{timestamp}.ini;
    # use section_name to record the output FASTA path in the config.
    saved_ini = encoder.save_config_file(section_name=out_fasta)
    out_ini = saved_ini  # returned path (timestamped)

    print(
        f"✓ Embedded {len(new_packets)} packet(s).  Pool saved to '{out_fasta}' / '{out_ini}'.",
        flush=True,
    )
    return 0


def _cmd_query(args: argparse.Namespace) -> int:
    """Search a FASTA pool for metadata sequences."""
    wanted = _resolve_metadata(
        getattr(args, "metadata", None), getattr(args, "metadata_file", None)
    )
    if not wanted:
        print("ERROR: no metadata sequences provided.", flush=True)
        return 1

    fasta = typing.cast(typing.Dict[str, str], load_fasta(args.fasta))
    total_seqs = len(fasta)

    results: typing.Dict[str, typing.List[str]] = {m: [] for m in wanted}
    for seq_id, dna in fasta.items():
        for meta in wanted:
            if meta in dna:
                results[meta].append(seq_id)

    found_any = False
    for meta, seq_ids in results.items():
        count = len(seq_ids)
        if count:
            found_any = True
        status = f"FOUND ({count} packet(s))" if count else "NOT FOUND"
        print(f"[{status}]  {meta}")
        if count and not args.count_only:
            for sid in seq_ids:
                print(f"        → {sid}")

    print(f"\nSearched {total_seqs} sequences in '{args.fasta}'.")
    return 0 if found_any else 2  # exit 2 = searched OK but not found


def _cmd_info(args: argparse.Namespace) -> int:
    """Show padding capacity for metadata embedding."""
    cfg_worker, decoder, semiautomatic_solver = _load_pool(args.ini)

    # Solve to populate GEPP
    decoder.solve()
    decoder.populate_header_chunk()

    tmp_gepp = semiautomatic_solver.decoder.GEPP.clone()
    tmp_gepp.solve()
    header_row = int(tmp_gepp.result_mapping[0][0])
    last_chunk_len_str = semiautomatic_solver.decoder.config_map.get("last_chunk_len_str", "I")
    header_chunk = HeaderChunk(
        Packet(
            tmp_gepp.b[header_row],
            {0},
            semiautomatic_solver.decoder.number_of_chunks,
            read_only=True,
        ),
        last_chunk_len_format=last_chunk_len_str,
        checksum_len_format=semiautomatic_solver.decoder.checksum_len_str,
    )

    offset_no_fn, length_no_fn = get_padding_area_offset(header_chunk, filename_is_free_area=False)
    offset_fn, length_fn = get_padding_area_offset(header_chunk, filename_is_free_area=True)

    chunk_size = int(semiautomatic_solver.decoder.config_map.get("chunk_size", len(tmp_gepp.b[0])))
    num_chunks = semiautomatic_solver.decoder.number_of_chunks
    num_packets = len(semiautomatic_solver.decoder.packets)
    filename = header_chunk.get_file_name()
    try:
        fn_str = filename.decode("utf-8", errors="replace")
    except AttributeError:
        fn_str = str(filename)

    print(f"Pool:              {args.ini}")
    print(f"Filename in pool:  {fn_str!r}  ({len(filename)} bytes)")
    print(f"Chunk size:        {chunk_size} bytes")
    print(f"Number of chunks:  {num_chunks}")
    print(f"Packets in pool:   {num_packets}")
    print()
    print("── Padding capacity (filename area kept) ─────────────────────────")
    print(f"  Padding offset:  {offset_no_fn} bytes")
    print(
        f"  Padding length:  {length_no_fn} bytes  →  max {length_no_fn * 4} DNA bases per metadata string"
    )
    print()
    print("── Padding capacity (filename area used as free area) ────────────")
    print(f"  Padding offset:  {offset_fn} bytes")
    print(
        f"  Padding length:  {length_fn} bytes  →  max {length_fn * 4} DNA bases per metadata string"
    )
    print()

    # Count packets covering chunk 0 (header) and last chunk
    n_header = sum(
        1
        for p in semiautomatic_solver.decoder.packets
        if not isinstance(p, str) and semiautomatic_solver.decoder.removeAndXorAuxPackets(p)[0]
    )
    n_last = sum(
        1
        for p in semiautomatic_solver.decoder.packets
        if not isinstance(p, str) and semiautomatic_solver.decoder.removeAndXorAuxPackets(p)[-1]
    )
    print(f"Header-chunk carrier candidates:    {n_header}")
    print(f"Last-chunk carrier candidates:      {n_last}")
    return 0


if __name__ == "__main__":
    import sys

    parsed_args = init_args()
    if parsed_args.command == "embed":
        sys.exit(_cmd_embed(parsed_args))
    elif parsed_args.command == "query":
        sys.exit(_cmd_query(parsed_args))
    elif parsed_args.command == "info":
        sys.exit(_cmd_info(parsed_args))
