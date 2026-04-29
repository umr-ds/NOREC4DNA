"""Common interface for degree distributions used by the encoders and decoders."""

import typing
from abc import ABC, abstractmethod

import numpy as np


class RandomWithSeed(typing.Protocol):
    def seed(self, seed: int) -> None: ...


class RandomWithChoice(RandomWithSeed, typing.Protocol):
    def choice(
        self,
        a: typing.Any,
        size: typing.Optional[int] = None,
        replace: bool = True,
        p: typing.Optional[typing.Sequence[float]] = None,
    ) -> typing.Any: ...


class Distribution(ABC):
    """Abstract base class for all supported packet-degree distributions."""

    def __init__(self) -> None:
        self.rng: RandomWithChoice = np.random
        self.pre_comp_dist: typing.List[float] = []
        self.S: typing.Optional[int] = None

    def get_config_string(self) -> str:
        return "Interface"

    @staticmethod
    def normalize(dist: typing.List[float]) -> typing.List[float]:
        return [float(i) / sum(dist) for i in dist]

    def get_distribution(self) -> typing.List[float]:
        return self.pre_comp_dist

    def get_size(self) -> typing.Optional[int]:
        return self.S

    def set_seed(self, seed: int) -> None:
        self.rng.seed(seed)

    @abstractmethod
    def update_number_of_chunks(self, num_chunks: int) -> None:
        pass  # implemented in subclasses

    @abstractmethod
    def getNumber(self, *args: typing.Any, **kwargs: typing.Any) -> int:
        pass
