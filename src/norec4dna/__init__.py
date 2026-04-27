"""Public package interface for NOREC4DNA."""

from importlib.metadata import PackageNotFoundError, version

__all__ = [
    "DecodePacket",
    "Decoder",
    "GEPP",
    "helper",
    "rules",
    "distributions",
    "LTEncoder",
    "RU10Encoder",
    "RU10Packet",
    "Encoder",
    "HeaderChunk",
    "LTBPDecoder",
    "LTDecoder",
    "OnlineBPDecoder",
    "OnlineDecoder",
    "OnlineEncoder",
    "OnlineAuxPacket",
    "OnlinePacket",
    "Packet",
    "ReedSolomonEncoder",
    "ReedSolomonDecoder",
    "ReedSolomonSuite",
    "RU10Decoder",
    "RU10BPDecoder",
    "RU10IntermediatePacket",
    "ErrorCorrection",
    "nocode",
    "reed_solomon_encode",
    "reed_solomon_decode",
    "crc32",
]

from . import ErrorCorrection
from . import ReedSolomonSuite
from . import distributions
from . import helper
from . import rules
from .DecodePacket import DecodePacket
from .Decoder import Decoder
from .Encoder import Encoder
from .ErrorCorrection import crc32, nocode, reed_solomon_decode, reed_solomon_encode
from .GEPP import GEPP
from .HeaderChunk import HeaderChunk
from .LTBPDecoder import LTBPDecoder
from .LTDecoder import LTDecoder
from .LTEncoder import LTEncoder
from .OnlineAuxPacket import OnlineAuxPacket
from .OnlineBPDecoder import OnlineBPDecoder
from .OnlineDecoder import OnlineDecoder
from .OnlineEncoder import OnlineEncoder
from .OnlinePacket import OnlinePacket
from .Packet import Packet
from .ReedSolomonSuite import ReedSolomonDecoder, ReedSolomonEncoder
from .RU10BPDecoder import RU10BPDecoder
from .RU10Decoder import RU10Decoder
from .RU10Encoder import RU10Encoder
from .RU10IntermediatePacket import RU10IntermediatePacket
from .RU10Packet import RU10Packet

__title__ = "norec4dna"
try:
    __version__ = version("norec4dna")
except PackageNotFoundError:
    __version__ = "0.2.0"
__author__ = "Michael Schwarz"
__copyright__ = "Copyright 2018 Michael Schwarz"
