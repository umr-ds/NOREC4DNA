# From MestreLion @ https://stackoverflow.com/a/38317060/9742086
from threading import Timer
from typing import Any, Callable, Optional


class RepeatedTimer(object):
    def __init__(self, interval: float, function: Callable[..., Any], *args: Any, **kwargs: Any):
        self._timer: Optional[Timer] = None
        self.interval = interval
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self.is_running = False
        self.start()

    def _run(self) -> None:
        self.is_running = False
        self.start()
        self.function(*self.args, **self.kwargs)

    def start(self) -> None:
        if not self.is_running:
            self._timer = Timer(self.interval, self._run)
            self._timer.start()
            self.is_running = True

    def stop(self) -> None:
        if self._timer is not None:
            self._timer.cancel()
        self.is_running = False
