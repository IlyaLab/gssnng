import logging

logger = logging.getLogger('gssnng')
logger.setLevel(logging.INFO)
logger.propagate = False
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %('
                              'messages)s')
handler.setFormatter(formatter)

if not logger.handlers:
    logger.addHandler(handler)
