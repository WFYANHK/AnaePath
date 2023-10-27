import logging

## setup logger format
logging.basicConfig(
    level="INFO",
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S")

logging.addLevelName(logging.WARNING,
                     f'\033[1m\x1b[33;20m{logging.getLevelName(logging.WARNING)}\033[1;0m')
logging.addLevelName(logging.CRITICAL,
                     f'\033[1m\x1b[31;20m{logging.getLevelName(logging.CRITICAL)}\033[1;0m')

logger = logging.getLogger(__name__)
