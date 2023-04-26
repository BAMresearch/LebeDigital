from lebedigital.openbis.interbis import Interbis
from pydantic import BaseModel, validator
from typing import Optional, Union
from pybis.entity_type import SampleType

class SampleBase(BaseModel):
    type: Union[SampleType, str] = ""
    space: str = ""
    project: str = ""
    experiment: str = ""
    parents: Optional[list[str]] = [] 
    children: Optional[list[str]] = [] 


