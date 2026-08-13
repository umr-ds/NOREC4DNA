import os
from collections import deque
from typing import IO, Deque, Dict, List, Optional, Set, Union

from .Decoder import Decoder
from .distributions import Distribution
from .ErrorCorrection import ErrorCorrectionCallable, nocode
from .HeaderChunk import HeaderChunk
from .helper.RU10Helper import from_true_false_list
from .OnlineAuxPacket import OnlineAuxPacket
from .OnlinePacket import OnlinePacket
from .Packet import Packet
from .RU10IntermediatePacket import RU10IntermediatePacket
from .RU10Packet import RU10Packet


class BPDecoder(Decoder):
    def __init__(
        self,
        file: Optional[str] = None,
        error_correction: ErrorCorrectionCallable = nocode,
        use_headerchunk: bool = True,
        static_number_of_chunks: Optional[int] = None,
        use_method: bool = False,
    ):
        super().__init__(file=file)
        self.debug = False
        self.isPseudo: bool = False
        self.file: Optional[str] = file
        self.use_method: bool = use_method
        self.f: Optional[IO[bytes]] = None
        if file is not None:
            self.isFolder = os.path.isdir(file)
            if not self.isFolder:
                try:
                    self.f = open(file, "rb")
                except Exception:
                    self.f = None
        self.correct: int = 0
        self.corrupt: int = 0
        self.number_of_chunks: int = 1000000
        self.headerChunk: Optional[HeaderChunk] = None
        self.decodedPackets: Set[Packet] = set()
        self.queue: Deque[Packet] = deque()
        self.symbol_to_packets: Dict[int, Set[Packet]] = {}
        self.solved_symbols: Dict[int, Packet] = {}
        self.pseudoCount: int = 0
        self.repairBlockNumbers: Dict[int, Set[int]] = {}
        self.s: int = -1
        self.h: int = -1
        self.numberOfDecodedAuxBlocks: int = 0
        self.dist: Optional[Distribution] = None
        self.EOF: bool = False
        self.counter: Dict[int, int] = {}
        self.count: bool = False
        self.error_correction: ErrorCorrectionCallable = error_correction
        self.use_headerchunk: bool = use_headerchunk
        self.static_number_of_chunks: Optional[int] = static_number_of_chunks
        self.auxBlocks: Dict[int, Union[RU10IntermediatePacket, OnlineAuxPacket]] = {}
        self.degreeToPacket: Dict[int, Set[Packet]] = {}
        self.ldpcANDhalf: Dict[int, RU10IntermediatePacket] = {}

    def addPacket(self, packet: Union[Packet, RU10Packet, OnlinePacket]) -> None:
        removed = self.removeAndXorAuxPackets(packet)
        packet.set_used_packets(set(from_true_false_list(removed)))
        if self.static_number_of_chunks is None:
            self.number_of_chunks = packet.get_total_number_of_chunks()

    def updatePackets(self, packet: Packet) -> bool:
        deg = packet.get_degree()
        if deg == 0:
            return self.is_decoded()

        # Un-XOR already solved symbols from this incoming packet first
        for sym_id in list(packet.get_used_packets()):
            if sym_id in self.solved_symbols:
                packet.xor_and_remove_packet(self.solved_symbols[sym_id])

        deg = packet.get_degree()
        if deg == 0:
            return self.is_decoded()

        if deg == 1:
            self.queue.append(packet)
        else:
            for sym_id in packet.get_used_packets():
                if sym_id not in self.symbol_to_packets:
                    self.symbol_to_packets[sym_id] = set()
                self.symbol_to_packets[sym_id].add(packet)

        # Peeling only: the expensive inactivation step (if any) is deferred to
        # the explicit solve() call so the per-packet path stays near-linear.
        return self._peel()

    def _peel(self) -> bool:
        """Run one belief-propagation peeling pass (near-linear).

        Solves every symbol that can be resolved from degree-1 packets.  This is
        the cheap, incremental solver used on the per-packet hot path; the more
        expensive inactivation step (if any) is deferred to ``solve()``.
        """
        while len(self.queue) > 0:
            deg1_pkt = self.queue.popleft()
            if deg1_pkt.get_degree() != 1:
                # Reduce by known solved symbols
                for sym_id in list(deg1_pkt.get_used_packets()):
                    if sym_id in self.solved_symbols:
                        deg1_pkt.xor_and_remove_packet(self.solved_symbols[sym_id])
                if deg1_pkt.get_degree() != 1:
                    if deg1_pkt.get_degree() > 1:
                        for sym_id in deg1_pkt.get_used_packets():
                            if sym_id not in self.symbol_to_packets:
                                self.symbol_to_packets[sym_id] = set()
                            self.symbol_to_packets[sym_id].add(deg1_pkt)
                    continue

            [sym_id] = deg1_pkt.get_used_packets()
            if sym_id >= self.number_of_chunks:
                continue

            if sym_id in self.solved_symbols:
                continue

            # Store solved symbol
            self.solved_symbols[sym_id] = deg1_pkt
            if sym_id == 0 and self.use_headerchunk and self.headerChunk is None:
                try:
                    last_chunk_fmt = "I"
                    if hasattr(self, "config_map") and self.config_map is not None:
                        last_chunk_fmt = str(self.config_map.get("last_chunk_len_str", "I"))
                    checksum_fmt = getattr(self, "checksum_len_str", "") or ""
                    self.headerChunk = HeaderChunk(
                        deg1_pkt,
                        last_chunk_len_format=last_chunk_fmt,
                        checksum_len_format=checksum_fmt,
                    )
                except Exception:
                    pass
            if sym_id != 0 or not self.use_headerchunk:
                self.decodedPackets.add(deg1_pkt)

            # Propagate symbol_id to all unreduced packets containing symbol_id
            affected_packets = self.symbol_to_packets.pop(sym_id, set())
            for other_pkt in affected_packets:
                if other_pkt is deg1_pkt:
                    # Never XOR a packet with itself: this would zero the just-solved symbol.
                    continue
                if sym_id in other_pkt.get_used_packets():
                    other_pkt.xor_and_remove_packet(deg1_pkt)
                    new_deg = other_pkt.get_degree()
                    if new_deg == 1:
                        self.queue.append(other_pkt)

        return self.is_decoded()

    def solve(self) -> bool:
        return self._peel()

    def removeAndXorAuxPackets(
        self, packet: Union[Packet, RU10Packet, OnlinePacket]
    ) -> Union[List[bool], Set[int]]:
        return []

    def is_decoded(self) -> bool:
        return len(self.decodedPackets) + (1 if self.headerChunk is not None else 0) >= self.number_of_chunks
