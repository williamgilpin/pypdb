"""Utility functions for requesting URLs over HTTP"""

from typing import Optional

import time
import requests
import warnings


def request_limited(url: str,
                    rtype: str = "GET",
                    num_attempts: int = 3,
                    sleep_time=0.5,
                    **kwargs) -> Optional[requests.models.Response]:
    """
    HTML request with rate-limiting base on response code


    Parameters
    ----------
    url : str
        The url for the request
    rtype : str
        The request type (oneof ["GET", "POST"])
    num_attempts : int
        In case of a failed retrieval, the number of attempts to try again
    sleep_time : int
        The amount of time to wait between requests, in case of
        API rate limits
    **kwargs : dict
        The keyword arguments to pass to the request

    Returns
    -------

    response : requests.models.Response
        The server response object. Only returned if request was successful,
        otherwise returns None.

    """

    if rtype not in ["GET", "POST"]:
        warnings.warn("Request type not recognized")
        return None

    total_attempts = 0
    while (total_attempts <= num_attempts):
        if rtype == "GET":
            response = requests.get(url, **kwargs)
        elif rtype == "POST":
            response = requests.post(url, **kwargs)

        if response.status_code == 200:
            return response

        if response.status_code == 429:
            curr_sleep = (1 + total_attempts) * sleep_time
            warnings.warn("Too many requests, waiting " + str(curr_sleep) +
                          " s")
            time.sleep(curr_sleep)
        elif 500 <= response.status_code < 600:
            warnings.warn("Server error encountered. Retrying")
        total_attempts += 1

    warnings.warn("Too many failures on requests. Exiting...")
    return None
