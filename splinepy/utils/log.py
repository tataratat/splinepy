"""splinepy/splinepy/log.py
Thin logging wrapper.
"""

import functools
import logging


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
    logger = logging.getLogger("splinepy")

    # level
    level = logging.DEBUG if debug else logging.INFO
    logger.setLevel(level)

    # format
    formatter = logging.Formatter(
        fmt="%(asctime)-15s %(name)s [%(levelname)s] %(message)s"
    )

    # apply format using stream handler
    # let's use only one stream handler so that calling configure multiple
    # times won't duplicate printing.
    new_handlers = list()
    for i, h in enumerate(logger.handlers):
        # we skip all the stream handler.
        if isinstance(h, logging.StreamHandler):
            continue

        # blindly keep other ones
        else:
            new_handlers.append(h)

    # add new stream handler
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(formatter)
    new_handlers.append(stream_handler)

    logger.handlers = new_handlers

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
