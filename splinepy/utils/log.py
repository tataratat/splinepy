"""splinepy/splinepy/log.py
Thin logging wrapper.
"""

import functools as _functools
import logging as _logging


def configure(debug=False, logfile=None):
    """
    Logging configurator. Can help you to set debug or info.
    Calling this will keep only one logging.StreamHandler.

    Parameters
    -----------
    debug: bool
    logfile: str

    Returns
    --------
    None
    """
    # logger
    logger = _logging.getLogger("splinepy")

    # level
    level = _logging.DEBUG if debug else _logging.INFO
    logger.setLevel(level)

    # format
    formatter = _logging.Formatter(
        fmt="%(asctime)-15s %(name)s [%(levelname)s] %(message)s"
    )

    # apply format using stream handler
    # let's use only one stream handler so that calling configure multiple
    # times won't duplicate printing.
    new_handlers = []
    for _i, h in enumerate(logger.handlers):
        # we skip all the stream handler.
        if isinstance(h, _logging.StreamHandler):
            continue

        # blindly keep other ones
        else:
            new_handlers.append(h)

    # add new stream handler
    stream_handler = _logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(formatter)
    new_handlers.append(stream_handler)

    logger.handlers = new_handlers

    # output logs
    if logfile is not None:
        file_logger_handler = _logging.FileHandler(logfile)
        logger.addHandler(file_logger_handler)


def debug(*log):
    """
    Debug logger.

    Parameters
    -----------
    *log: str

    Returns
    --------
    None
    """
    logger = _logging.getLogger("splinepy")
    logger.debug(" ".join(map(str, log)))


def info(*log):
    """
    Info logger.

    Parameters
    -----------
    *log: str

    Returns
    --------
    None
    """
    logger = _logging.getLogger("splinepy")
    logger.info(" ".join(map(str, log)))


def warning(*log):
    """
    warning logger.

    Parameters
    -----------
    *log: str

    Returns
    --------
    None
    """
    logger = _logging.getLogger("splinepy")
    logger.warning(" ".join(map(str, log)))


def prepend_log(message, log_func):
    """
    Prepend message before a logging function.

    Parameters
    ----------
    message: str
    log_func: function
      one of the following - {info, debug, warning}

    Returns
    -------
    prepended: function
    """
    return _functools.partial(log_func, message)
