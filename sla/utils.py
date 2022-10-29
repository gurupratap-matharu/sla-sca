import json
import logging

logger = logging.getLogger("utils")


def write_to_local(response, filename):
    """
    Writes a GEE object to the local filesystem in JSON format.

    Be cautious as this method will call the getInfo() method on the GEE
    object to retrieve all the results from the EE server.
    So it can be computationally expensive and block all further execution.
    """

    with open(filename, "w", encoding="utf-8") as f:
        f.write(json.dumps(response.getInfo()))

    logger.info(f"Results written to {filename}")
