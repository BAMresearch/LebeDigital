import logging
from math import isnan
from sys import exit
from typing import Iterable

import pandas as pd
import yaml
from termcolor import colored

from lebedigital.openbis.interbis import Interbis
from pybis.sample import Sample


class ISample(Sample):
    pass
