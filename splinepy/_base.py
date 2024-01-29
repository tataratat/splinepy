"""splinepy/splinepy/_base.py

Base type of splinepy
"""

from splinepy.utils import log as _log


class SplinepyBase:
    __slots__ = ()

    def __init_subclass__(cls, *args, **kwargs):
        """
        Add logger shortcut.
        """
        super().__init_subclass__(*args, **kwargs)
        cls._logi = _log.prepend_log("<" + cls.__qualname__ + ">", _log.info)
        cls._logd = _log.prepend_log("<" + cls.__qualname__ + ">", _log.debug)
        cls._logw = _log.prepend_log(
            "<" + cls.__qualname__ + ">", _log.warning
        )
        return super().__new__(cls)
