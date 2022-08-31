"""splinepy/splinepy/log.py
Thin logging wrapper.
"""

import logging


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
    logger = logging.getLogger()
    if debug:
        logger.setLevel(logging.DEBUG)

    else:
        logger.setLevel(logging.INFO)

    if logfile is not None:
        file_logger_handler = logging.FileHandler(logfile)


def debug(*log):
    """
    Debug logger.
    Parameters
    -----------
    *log: *str
    Returns
    --------
    None
    """
    logging.debug(" ".join(map(str, log)))


def info(*log):
    """
    Info logger.
    Parameters
    -----------
    *log: *str
    Returns
    --------
    None
    """
    logging.info(" ".join(map(str, log)))


def warning(*log):
    """
    warning logger.
    Parameters
    -----------
    *log: *str
    Returns
    --------
    None
    """
    logging.warning(" ".join(map(str, log)))