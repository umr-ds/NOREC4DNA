import argparse
import io
import logging
import random
import struct
import typing

import numpy
import numpy as np
from python_utils.types import deprecated

from NOREC4DNA.norec4dna import RU10Decoder
from norec4dna.helper.quaternary2Bin import tranlate_quat_to_byte

from norec4dna.helper.bin2Quaternary import string2QUATS

from NOREC4DNA.ConfigWorker import ConfigReadAndExecute
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.RU10Packet import RU10Packet
from .invivo_window_decoder import load_fasta
from norec4dna import RU10Encoder
from norec4dna.Packet import Packet
from norec4dna.helper.helper import calc_crc

from semi_automatic_reconstruction_toolkit import SemiAutomaticReconstructionToolkit
import imagehash

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def get_padding_area_offset(header_chunk: HeaderChunk, filename_is_free_area=False) -> typing.Tuple[int, int]:
    """
    Returns the position and length of the unused padding area from the header chunk
    @header_chunk: the header chunk to extract from
    @filename_is_free_area:
    @returns: Tuple containing the start offset of the padding area and the padding_length
    """
    padding_offset = struct.calcsize(header_chunk.last_chunk_len_format + header_chunk.last_chunk_len_format)
    if filename_is_free_area:
        padding_offset += len(
            header_chunk.get_file_name()) - 1  # -1 as we must ensure 1 byte contains 0x00 (zero terminated string!)
    padding_length = len(header_chunk.data) - padding_offset
    return padding_offset, padding_length


# Header should contain:
#   - # changed chunks
#   - hash of new file
#   - reference (e.g. by seed) to an artificial packet containing all packets (by seed) that were changed -> these would be inserted FIRST during decoding! (optional) -> this concept could also be used for inserting DIFF information (increasing and decreasing the size)

def update_header(encoder, semiautomatic_solver: SemiAutomaticReconstructionToolkit,
                  new_hash: None, new_filename: None, new_version: None):
    # new_header = HeaderChunk(packet)
    if new_filename is None:
        new_filename = ""

    old_header = HeaderChunk(
        Packet(semiautomatic_solver.decoder.GEPP.b[semiautomatic_solver.decoder.GEPP.result_mapping[0]], {0},
               semiautomatic_solver.decoder.number_of_chunks, read_only=True),
        last_chunk_len_format="I", checksum_len_format=semiautomatic_solver.decoder.checksum_len_str)
    new_header: HeaderChunk = old_header  # TODO: create new header with updated values
    tmp = np.frombuffer(encoder.chunks[0][1:], dtype=np.uint8).reshape(-1)
    checksum = calc_crc(io.BytesIO(
        tmp[:len(tmp) - old_header.last_chunk_length]))  # TODO!
    additional_payload = ""  # TODO: add metadata or special information and pad it to fill the chunk.
    old_header.update_header(new_filename, checksum, additional_payload)
    encoder.chunks.insert(0, np.frombuffer(new_header.data, dtype=np.uint8))
    pass


def calculate_new_last_chunk_packet(packet: RU10Packet, last_chunk_length: int, wanted_metadata: bytes) -> RU10Packet:
    assert len(wanted_metadata) <= len(packet.data) - last_chunk_length
    # wanted_metadata = np.zeros_like(np_header_chunk)
    # metadata_array = np.zeros_like(np_header_chunk)
    start_idx = len(packet.data) - len(wanted_metadata)
    # metadata_array[start_idx:] = np.frombuffer(wanted_metadata, dtype=np.bool)
    np_packet = numpy.frombuffer(packet.data, dtype=np.uint8).copy()
    logger.error(np_packet)
    logger.error("".join(string2QUATS(np_packet.tobytes())))
    np_packet[start_idx:] = np.frombuffer(wanted_metadata, dtype=np.uint8)
    logger.error(np_packet)
    logger.error("".join(string2QUATS(np_packet.tobytes())))

    out_packet = RU10Packet(np_packet, used_packets=packet.used_packets,
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
                            id_spacing=packet.id_spacing)
    return out_packet


def calculate_new_header_packet(packet: RU10Packet, header_chunk: HeaderChunk, wanted_metadata: bytes) -> RU10Packet:
    padding_offset, padding_length = get_padding_area_offset(header_chunk, filename_is_free_area=True)
    assert len(wanted_metadata) <= padding_length
    np_header_chunk = numpy.frombuffer(header_chunk.data, dtype=np.uint8)
    # wanted_metadata = np.zeros_like(np_header_chunk)
    # metadata_array = np.zeros_like(np_header_chunk)
    start_idx = len(np_header_chunk) - len(wanted_metadata)
    # metadata_array[start_idx:] = np.frombuffer(wanted_metadata, dtype=np.bool)
    np_packet = numpy.frombuffer(packet.data, dtype=np.uint8).copy()
    logger.error(np_packet)
    logger.error("".join(string2QUATS(np_packet.tobytes())))
    np_packet[start_idx:] = np.frombuffer(wanted_metadata, dtype=np.uint8)
    logger.error(np_packet)
    logger.error("".join(string2QUATS(np_packet.tobytes())))

    out_packet = RU10Packet(np_packet, used_packets=packet.used_packets,
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
                            id_spacing=packet.id_spacing)
    return out_packet


def create_new_metadata_packets(semiautomatic_solver: SemiAutomaticReconstructionToolkit, metadata: list[str],
                                unwanted_metadata: list[str], replace_packets_with_header: bool = True,
                                replace_packets_with_last_chunk: bool = False):
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
    if not (replace_packets_with_header and replace_packets_with_last_chunk):
        raise RuntimeError(
            "At least one of replace_packets_with_header or replace_packets_with_last_chunk must be set!")
    new_packets = []
    packets_using_header = []
    packets_using_last_chunk = []
    if len(semiautomatic_solver.decoder.packets) == 0:
        raise ValueError("Decoder does not contain any packets!")
    tmp_gepp = semiautomatic_solver.decoder.GEPP.clone()
    tmp_gepp.solve()
    header_row = tmp_gepp.result_mapping[0]
    header_chunk = HeaderChunk(
        Packet(tmp_gepp.b[header_row], {0}, semiautomatic_solver.decoder.number_of_chunks, read_only=True),
        last_chunk_len_format=semiautomatic_solver.last_chunk_len_format,
        checksum_len_format=semiautomatic_solver.decoder.checksum_len_str)
    added_metadata_str: set = set()
    for packet in semiautomatic_solver.decoder.packets:
        if type(packet) is str:
            # Raw DNA sequence as checksum / RS failed - indicator for metadata influence.. -> UNLESS we update the Checksum as well
            print("Invalid packet...")
            dna_str = packet
        else:  # RU10Packet or any other type:
            dna_str = string2QUATS(packet.data)
        for metadata_str in metadata:
            if metadata_str in dna_str:
                logger.warn(f"Found wanted metadata string {metadata_str} in a packet")

        for unwanted_metadata_str in unwanted_metadata:
            if unwanted_metadata_str in dna_str:
                logger.error(f"Found unwanted metadata string {unwanted_metadata_str} in a packet! "
                             f"This will lead to false-positives during querying!")
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
            new_packet = calculate_new_header_packet(packet, header_chunk, tranlate_quat_to_byte(metadata_str))
        else:
            new_packet = calculate_new_last_chunk_packet(packet, header_chunk.last_chunk_length,
                                                         tranlate_quat_to_byte(metadata_str))
        new_packets.append(new_packet)
        added_metadata_str.add(metadata_str)

    return new_packets

@deprecated("Use encoder_from_decoder instead!")
def encoder_from_sart(semiautomatic_solver: SemiAutomaticReconstructionToolkit,
                      config_worker: ConfigReadAndExecute, rules=None) -> RU10Encoder:
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
    ru10_encoder = RU10Encoder(semiautomatic_solver.decoder.headerChunk.file_name.decode(), number_of_chunks,
                               dist, semiautomatic_solver.decoder.use_headerchunk,
                               None, 0, rules,
                               semiautomatic_solver.decoder.error_correction, "I", crc_len_format,
                               number_of_chunks_len_format, id_len_format,
                               save_number_of_chunks_in_packet, False, "", "",
                               1.0, True, "", xor_by_seed,
                               mask_id, id_spacing, last_chunk_len_str)
    ru10_encoder.chunk_size = chunk_size  # avoid overwriting the number of chunks by postponing the chunk_size setup
    ru10_encoder.checksum_len_str = checksum_len_str
    ru10_encoder.checksum = semiautomatic_solver.decoder.headerChunk.checksum
    # remove the header and any decoded rows after the last chunk. Padding will be auto-applied as we use GEPP.b
    ru10_encoder.chunks = [x for x in semiautomatic_solver.decoder.GEPP.b[
        0:number_of_chunks]]  # TODO: check for off-by-one due to header
    ru10_encoder.generate_intermediate_blocks()  # missing command from the prepare block
    # (we MUST not call prepare() as it would add an additional header!)
    for packet in semiautomatic_solver.decoder.packets:
        packet.packed_used_packets = packet.prepare_and_pack()
        packet.packed = packet.calculate_packed_data()
    ru10_encoder.encodedPackets = set(semiautomatic_solver.decoder.packets.copy())

    # ru10_encoder.chunks.insert(0,header)
    # ru10_encoder.generate_intermediate_blocks()

    return ru10_encoder

def encoder_from_decoder(decoder: RU10Decoder,
                      config_worker: ConfigReadAndExecute, rules=None) -> RU10Encoder:
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
    ru10_encoder = RU10Encoder(decoder.headerChunk.file_name.decode(), number_of_chunks,
                               dist, decoder.use_headerchunk,
                               None, 0, rules,
                               decoder.error_correction, "I", crc_len_format,
                               number_of_chunks_len_format, id_len_format,
                               save_number_of_chunks_in_packet, False, "", "",
                               1.0, True, "", xor_by_seed,
                               mask_id, id_spacing, last_chunk_len_str)
    ru10_encoder.chunk_size = chunk_size  # avoid overwriting the number of chunks by postponing the chunk_size setup
    ru10_encoder.checksum_len_str = checksum_len_str
    ru10_encoder.checksum = decoder.headerChunk.checksum
    # remove the header and any decoded rows after the last chunk. Padding will be auto-applied as we use GEPP.b
    ru10_encoder.chunks = [x for x in decoder.GEPP.b[
        0:number_of_chunks]]  # TODO: check for off-by-one due to header
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
    output = []
    fasta = load_fasta(file)
    for name, meta in fasta.items():
        name, amount = name.split("_")
        output += int(amount) * [meta]
    return output


def replace_packets(encoder: RU10Encoder, new_packet_list: typing.List[RU10Packet], keep_percentage: float = 0.0,
                    keep_one_header=False, keep_one_last_chunk=False):
    """
    replace the matching packets (same seed)
    @encoder: encoder instance to update
    @new_packet_list: list of new packets to insert
    @keep_percentage: percentage of packets to
    """
    # TODO: we _MUST_ make sure that SOME of the original packets with the header are kept - otherwise the filename gets lost!
    # TODO: we _MUST_ make sure that each selected packet with the header can be used to decode WITHOUT any other modified packets - otherwise we may end up with a header that is constructed using two modified packets making it impossible or very hard to resolve the deltas
    # - or: calculate all deltas as we did for the automatic reconstruction to detect the unique error deltas and the possible packets to find the delta of each packet!
    num_packets_to_keep: int = int(np.ceil(len(new_packet_list) * keep_percentage))
    new_packets_seeds = [x.id for x in new_packet_list]
    packets_to_remove: typing.List[RU10Packet] = []
    # find all packet with same seed as the new packets
    for packet in encoder.encodedPackets:
        if packet.id in new_packets_seeds:
            packets_to_remove.append(packet)
    random.shuffle(packets_to_remove)
    # remove "unwanted" packets from encoder:
    # if len(packets_to_remove) > num_packets_to_keep:
    #    for packet in packets_to_remove[int(num_packets_to_keep):]:
    #        encoder.encodedPackets.remove(packet)
    num_packets_to_keep_offset = 0
    if keep_one_header:
        special_packet_to_keep = -1
        for packet_to_keep in range(num_packets_to_keep):
            if semiautomatic_solver.decoder.removeAndXorAuxPackets(packets_to_remove[packet_to_keep])[0]:
                special_packet_to_keep = packet_to_keep
                break
        for potential_packet_to_keep in range(num_packets_to_keep, len(packets_to_remove)):
            if semiautomatic_solver.decoder.removeAndXorAuxPackets(packets_to_remove[potential_packet_to_keep])[0]:
                special_packet_to_keep = potential_packet_to_keep
                break
        if special_packet_to_keep == -1:
            raise RuntimeError("No packet with header chunk found to keep!")
        # swap first with current packet_to_keep:
        packets_to_remove[0], packets_to_remove[special_packet_to_keep] = packets_to_remove[
            special_packet_to_keep], packets_to_remove[0]
        num_packets_to_keep_offset += 1

    if keep_one_last_chunk:
        special_packet_to_keep = -1
        for packet_to_keep in range(num_packets_to_keep):
            if semiautomatic_solver.decoder.removeAndXorAuxPackets(packets_to_remove[packet_to_keep])[-1]:
                special_packet_to_keep = packet_to_keep
                break
        for potential_packet_to_keep in range(num_packets_to_keep, len(packets_to_remove)):
            if semiautomatic_solver.decoder.removeAndXorAuxPackets(packets_to_remove[potential_packet_to_keep])[-1]:
                special_packet_to_keep = potential_packet_to_keep
                break
        if special_packet_to_keep == -1:
            raise RuntimeError("No packet with last chunk found to keep!")
        # swap second with current packet_to_keep:
        packets_to_remove[1], packets_to_remove[special_packet_to_keep] = packets_to_remove[
            special_packet_to_keep], packets_to_remove[1]
        num_packets_to_keep_offset += 1
    if num_packets_to_keep < num_packets_to_keep_offset: num_packets_to_keep = num_packets_to_keep_offset
    to_remove = set(packets_to_remove[int(num_packets_to_keep):])
    to_remove_ids = {id(p) for p in to_remove}
    # Rebuild without relying on equality lookups during remove
    encoder.encodedPackets = {p for p in encoder.encodedPackets if id(p) not in to_remove_ids}
    skipped_packets = 0
    # add all new packets:
    for packet in new_packet_list:
        if packet in encoder.encodedPackets:
            logger.error(f"Packet already in encoded Packets - this packet will be skipped: {packet}")
            skipped_packets += 1
        encoder.encodedPackets.add(packet)
    logger.warning(f"There were a total of {skipped_packets} packets containing duplicates. "
                   f"The the actual amount of copies may be different from what you requested. "
                   f"This usually indicates that there the requested number of copies is larger that the total number of packets containing the header chunk.")


def init_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ini", metavar="ini", type=str, help="config file (ini)",
                        default="/home/michael/Code/DR4DNA/eval/sleeping_beauty_no_error.ini")
    wanted_arg_group = parser.add_mutually_exclusive_group(required=True)
    wanted_arg_group.add_argument("--wanted_metadata_file", metavar="wantedmetafile", type=str,
                                  help="file containing unwanted metadata in the fasta format, description should contain a multiplier and the corresponding metadata ")
    wanted_arg_group.add_argument("--wanted_metadata", metavar="wantedmeta", type=str,
                                  help="comma-separated list of metadata DNA sequences")
    unwanted_arg_group = parser.add_mutually_exclusive_group(required=False)
    unwanted_arg_group.add_argument("--unwanted_metadata_file", metavar="unwantedmetafile", type=str,
                                    help="file containing unwanted metadata in the fasta format, description should contain a multiplier and the corresponding metadata ")
    unwanted_arg_group.add_argument("--unwanted_metadata", metavar="wantedmeta", type=str,
                                    help="comma-separated list of metadata DNA sequences")
    return parser.parse_args()


if __name__ == '__main__':
    # LES detektion der wiedersprüchlichen zeilen durch augmented matrix und diese mittels gauss lösen: [0 0 0 | x (x!=0)] -> beteiligte Zeilen sind verantwortlich!
    #
    parsed_args = init_args()

    ini_file = parsed_args.ini

    if parsed_args.wanted_metadata_file is not None:
        wanted_metadata = parse_metadata_file(parsed_args.wanted_metadata_file)
    else:
        wanted_metadata = parsed_args.wanted_metadata.split(",")

    if parsed_args.unwanted_metadata_file is not None:
        unwanted_metadata = parse_metadata_file(parsed_args.unwanted_metadata_file)
    else:
        unwanted_metadata = parsed_args.unwanted_metadata.split(",")

    cfg_worker = ConfigReadAndExecute(ini_file)
    x = cfg_worker.execute(return_decoder=True, skip_solve=True)[0]
    semiautomatic_solver = SemiAutomaticReconstructionToolkit(x)

    encoder = encoder_from_sart(semiautomatic_solver, cfg_worker)

    encoder.save_packets_fasta("1test_out.fasta", "", False)

    new_packets = create_new_metadata_packets(semiautomatic_solver, wanted_metadata, unwanted_metadata, True, True)

    replace_packets(encoder, new_packets, 0.0, True, True)

    logger.error(encoder.encodedPackets)

    encoder.save_packets_fasta("test_out.fasta", "", False)
    encoder.save_config_file(section_name="test_out.fasta")
    logger.warning(semiautomatic_solver.decoder.packets)
    # semiautomatic_solver.get_file_as_bytes()
