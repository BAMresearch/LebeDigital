import logging
import os
import sys
from getpass import getpass

import pandas as pd
from pybis import Openbis
from pybis.sample import Sample
import re

"""
This class will be very similar to pybis until I figure out what i can take from pybis and what I have to write myself
"""


class RequestsGenerator():
    """
    Here trying to hard code a sample fetch dict
    """

    def _create_get_request(self, method_name, entity, permids, options, foType):

        if not isinstance(permids, list):
            permids = [permids]

            type = f"as.dto.{entity.lower()}.id.{entity.capitalize()}"
            search_params = []
            for permid in permids:
                # decide if we got a permId or an identifier
                match = re.match("/", permid)
                if match:
                    search_params.append(
                        {"identifier": permid, "@type": type + "Identifier"}
                    )
                else:
                    search_params.append(
                        {"permId": permid, "@type": type + "PermId"})

            fo = {"@type": foType}
            for option in options:
                fo[option] = get_fetchoption_for_entity(option)

            request = {
                "method": method_name,
                "params": [self.token, search_params, fo],
            }
            return request

    """
    Method for fetching samples
    "sample": {
        "@type": "as.dto.sample.fetchoptions.SampleFetchOptions",
        "type": {"@type": "as.dto.sample.fetchoptions.SampleTypeFetchOptions"},
    }
    """
