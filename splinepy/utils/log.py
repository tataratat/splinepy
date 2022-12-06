"""splinepy/splinepy/log.py
Thin logging wrapper.
"""

import logging
import functools


def configure(debug=False, logfile=None):
    """
    Logging configurator.

    Parameters
    -----------
    debug: bool
    logfile: str

    Returns
    --------
    None
    """
    # logger
    logger = logging.getLogger("splinepy")

    if debug:
        logger.setLevel(logging.DEBUG)

    else:
        logger.setLevel(logging.INFO)

    # format
    formatter = logging.Formatter(
            fmt="%(asctime)-15s %(name)s [%(levelname)s] %(message)s"
    )

    # apply format using stream handler
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)

    # output logs
    if logfile is not None:
        file_logger_handler = logging.FileHandler(logfile)
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
    logger = logging.getLogger("splinepy")
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
    logger = logging.getLogger("splinepy")
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
    logger = logging.getLogger("splinepy")
    logger.warning(" ".join(map(str, log)))


def prepend_log(message, log_func):
    """
    Prepend message before a logging function.

    Parameters
    ----------
    messgae: str
    log_func: function
      one of the followings - {info, debug, warning}

    Returns
    -------
    prepended: function
    """
    return functools.partial(log_func, message)
