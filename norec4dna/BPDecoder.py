import os
from collections import deque
from typing import Any, Callable, Deque, Dict, List, Optional, Set, Tuple, Union

from norec4dna import Decoder, HeaderChunk
from norec4dna.distributions import Distribution
from norec4dna.ErrorCorrection import nocode
from norec4dna.helper.RU10Helper import from_true_false_list
from norec4dna.OnlineAuxPacket import OnlineAuxPacket
from norec4dna.OnlinePacket import OnlinePacket
from norec4dna.Packet import Packet
from norec4dna.RU10IntermediatePacket import RU10IntermediatePacket
from norec4dna.RU10Packet import RU10Packet


class BPDecoder(Decoder):
    def __init__(
        self,
        file: Optional[str] = None,
        error_correction: Callable[[bytes], bytes] = nocode,  # type: ignore[assignment]
        use_headerchunk: bool = True,
        static_number_of_chunks: Optional[int] = None,
        use_method: bool = False,
    ):
        super().__init__()
        self.debug = False
        self.isPseudo: bool = False
        self.file: Optional[str] = file
        self.use_method: bool = use_method
        if file is not None:
            self.isFolder = os.path.isdir(file)
            if not self.isFolder:
                self.f = open(self.file, "rb")
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
        self.error_correction: Callable[[bytes], bytes] = error_correction
        self.use_headerchunk: bool = use_headerchunk
        self.static_number_of_chunks: Optional[int] = static_number_of_chunks
        self.auxBlocks: Dict[int, Union[RU10IntermediatePacket, OnlineAuxPacket]] = {}
        self.degreeToPacket: Dict[int, Set[Packet]] = {}
        self.ldpcANDhalf: Dict[int, RU10IntermediatePacket] = {}

    def addPacket(self, packet: Union[Packet, RU10Packet, OnlinePacket]) -> None:
        removed = self.removeAndXorAuxPackets(packet)
        packet.set_used_packets(set(from_true_false_list(removed)))
        if (not packet.get_degree() in self.degreeToPacket) or (
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
        pass  # type: ignore[return-value]

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
        fin: bool = False
        lookup: List[int] = [i for i in self.degreeToPacket.keys() if packet.get_degree() < i]
        for i in lookup:
            if not isinstance(self.degreeToPacket[i], set):
                self.degreeToPacket[i] = set()
            for p in self.degreeToPacket[i].copy():
                p_used = p.get_used_packets()
                pack_used = packet.get_used_packets()
                if len(pack_used) < len(p_used) and pack_used.issubset(p_used):
                    self.degreeToPacket[i].remove(p)
                    degree = self.compareAndReduce(p, packet)
                    if isinstance(degree, bool) and degree is True:
                        return degree
        degree: int = packet.get_degree()
        lookup = [i for i in self.degreeToPacket.keys() if packet.get_degree() > i]
        for i in lookup:
            if not isinstance(self.degreeToPacket[i], set):
                self.degreeToPacket[i] = set()
            for p in self.degreeToPacket[i].copy():
                p_used = p.get_used_packets()
                pack_used = packet.get_used_packets()
                if len(pack_used) > len(p_used) and p_used.issubset(pack_used):
                    try:
                        self.degreeToPacket[degree].remove(packet)
                        degree = self.compareAndReduce(packet, p)
                        if isinstance(degree, bool) and degree is True:
                            return degree
                    except Exception:
                        continue
                    # If we already reduced a Packet with the same used_packets, there is no need to do it again
        return fin or self.is_decoded()
