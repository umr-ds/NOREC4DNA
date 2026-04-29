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
                self.f = open(file, "rb")
        self.correct: int = 0
        self.corrupt: int = 0
        self.number_of_chunks: int = 1000000
        self.headerChunk: Optional[HeaderChunk] = None
        self.decodedPackets: Set[Packet] = set()
        self.queue: Deque[Packet] = deque()
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
        if (packet.get_degree() not in self.degreeToPacket) or (
            not isinstance(self.degreeToPacket[packet.get_degree()], set)
        ):
            self.degreeToPacket[packet.get_degree()] = set()
        if self.static_number_of_chunks is None:
            self.number_of_chunks = packet.get_total_number_of_chunks()
        self.degreeToPacket[packet.get_degree()].add(packet)
        # Correct

    def updatePackets(self, packet: Packet) -> bool:
        if (
            packet.get_degree() == 1
            and next(iter(packet.get_used_packets())) < self.number_of_chunks
            and (next(iter(packet.get_used_packets())) != 0 or not self.use_headerchunk)
        ):
            # Directly add Packets that are degree == 1 (except for HeaderPacket)
            self.decodedPackets.add(packet)
            return self.is_decoded()
        self.queue.append(packet)
        return self.solve()

    def solve(self) -> bool:
        finished = False
        while len(self.queue) > 0 and not finished:
            finished = self.reduceAll(self.queue.popleft())
        return finished

    def removeAndXorAuxPackets(self, packet: Union[Packet, RU10Packet, OnlinePacket]) -> List[bool]:
        # Abstract method - implemented in subclasses
        # Return empty list as default (no packets removed)
        return []

    def compareAndReduce(self, packet: Packet, other: Packet) -> Union[bool, int]:
        if self.file is None:
            packet.remove_packets(other.get_used_packets())
        else:
            packet.xor_and_remove_packet(other)
        degree = packet.get_degree()
        if (degree not in self.degreeToPacket) or (
            not isinstance(self.degreeToPacket[degree], set)
        ):
            self.degreeToPacket[degree] = set()
        if degree == 1:
            [x] = packet.get_used_packets()  # Unpacking -> Fastest way to extract Element from Set
            if x > self.number_of_chunks:  # we got a new AUX-Packet
                raise RuntimeError("this should not have happened!")
            else:
                if x != 0 or not self.use_headerchunk:
                    self.decodedPackets.add(packet)  # Add Packet to decoded Packets
        self.degreeToPacket[degree].add(packet)
        if self.is_decoded():
            return True
        self.queue.append(packet)
        return degree

    def reduceAll(self, packet: Packet) -> bool:
        # lookup all packets for this to solve with ( when this packet has a subset of used Packets)
        if self._reduce_larger_degree_packets(packet):
            return True
        if self._reduce_smaller_degree_packets(packet):
            return True
        return self.is_decoded()

    def _normalize_degree_packets(self, degree: int):
        if not isinstance(self.degreeToPacket[degree], set):
            self.degreeToPacket[degree] = set()
        return self.degreeToPacket[degree]

    def _reduce_larger_degree_packets(self, packet: Packet) -> bool:
        lookup: List[int] = [i for i in self.degreeToPacket.keys() if packet.get_degree() < i]
        for degree in lookup:
            degree_packets = self._normalize_degree_packets(degree)
            for candidate in degree_packets.copy():
                packet_used = packet.get_used_packets()
                candidate_used = candidate.get_used_packets()
                if len(packet_used) < len(candidate_used) and packet_used.issubset(candidate_used):
                    degree_packets.remove(candidate)
                    reduced_degree = self.compareAndReduce(candidate, packet)
                    if isinstance(reduced_degree, bool) and reduced_degree is True:
                        return True
        return False

    def _reduce_smaller_degree_packets(self, packet: Packet) -> bool:
        degree: int = packet.get_degree()
        lookup = [i for i in self.degreeToPacket.keys() if packet.get_degree() > i]
        for lookup_degree in lookup:
            self._normalize_degree_packets(lookup_degree)
            for candidate in self.degreeToPacket[lookup_degree].copy():
                packet_used = packet.get_used_packets()
                candidate_used = candidate.get_used_packets()
                if len(packet_used) > len(candidate_used) and candidate_used.issubset(packet_used):
                    try:
                        self.degreeToPacket[degree].remove(packet)
                        degree = self.compareAndReduce(packet, candidate)
                        if isinstance(degree, bool) and degree is True:
                            return True
                    except Exception:
                        continue
        return False
