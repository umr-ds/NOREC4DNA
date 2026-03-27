import typing
from typing import Any, List, Optional

import progressbar
from norec4dna.Packet import Packet
from abc import ABC, abstractmethod

class Decoder(ABC):
    def __init__(self, file: Optional[str] = None):
        self.file: Optional[str] = file
        self.read_all_before_decode: bool = False
        self.isFolder: bool = False
        self.isZip: bool = False
        self.progress_bar: Optional[progressbar.ProgressBar] = None
        self.number_of_chunks: int = 1000000
        self.isPseudo: bool = False

    @staticmethod
    def create_progress_bar(max_value: int) -> progressbar.ProgressBar:
        widgets: List[Any] = [
            progressbar.Percentage(),
            progressbar.Bar(),
            " Correct: ",
            progressbar.Counter(),
            ", ",
            progressbar.Variable("Corrupt"),
            ", ",
            progressbar.AdaptiveETA(),
            " ",
            progressbar.Timer(),
        ]
        return progressbar.ProgressBar(
            max_value=max_value, widgets=widgets, max_error=False, redirect_stdout=True
        ).start()

    @classmethod
    def pseudo_decoder(
        cls,
        number_of_chunks: Optional[int] = None,
        read_all_before_decode: bool = False,
    ) -> "Decoder":

        def warn_pseudo() -> None:
            print("This method is not allowed while using pseudo Decoder")

        pseudo = cls(None)
        pseudo.decodeFile = warn_pseudo  # type: ignore[assignment]
        pseudo.decode = warn_pseudo  # type: ignore[assignment]
        pseudo.getNextValidPacket = warn_pseudo  # type: ignore[assignment]
        pseudo.saveDecodedFile = warn_pseudo  # type: ignore[assignment]
        pseudo.read_all_before_decode = read_all_before_decode
        if number_of_chunks is not None:
            pseudo.number_of_chunks = number_of_chunks
        pseudo.isPseudo = True
        return pseudo
    
    @abstractmethod
    def input_new_packet(self, packet: Packet) -> bool:
        pass  # implemented in subclasses

    @abstractmethod
    def is_decoded(self) -> bool:
        pass  # implemented in subclasses

    @abstractmethod
    def solve(self) -> bool:
        pass  # implemented in subclasses

    @abstractmethod
    def decodeFolder(self, *args: Any, **kwargs: Any) -> Optional[Any]:
        pass  # implemented in subclasses

    @abstractmethod
    def decodeFile(self, *args: Any, **kwargs: Any) -> Optional[Any]:
        pass  # implemented in subclasses

    @abstractmethod
    def decodeZip(self, *args: Any, **kwargs: Any) -> Optional[Any]:
        pass  # implemented in subclasses

    def set_read_all_before_decode(self, do: bool) -> None:
        self.read_all_before_decode = do

    def decode(self, *args: Any, **kwargs: Any) -> Optional[Any]:
        if self.isFolder:
            return self.decodeFolder(*args, **kwargs)
        elif self.isZip:
            return self.decodeZip(*args, **kwargs)
        else:
            return self.decodeFile(*args, **kwargs)
