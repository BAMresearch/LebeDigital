import functools
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()


class UploadLogger:
    def __init__(self):
        logging.basicConfig(level=logging.DEBUG)

    def get_logger(self, name=None):
        return logging.getLogger(name)


def log(func=None):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            result = func(*args, **kwargs)
            return result
        except Exception as e:
            logger.exception(f"Exception raised in {func.__name__}.\nException: {str(e)}")
    return wrapper
