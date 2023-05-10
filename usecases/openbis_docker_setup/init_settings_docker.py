from lebedigital.openbis.interbis import Interbis
from pathlib import Path
import json
import argparse

settings_path = Path("./init_settings.json")

argParser = argparse.ArgumentParser()
argParser.add_argument("-l", "--login", help="openbis login")
argParser.add_argument("-p", "--password", help="openbis password")
argParser.add_argument("-u", "--url", help="openbis url address")

args = argParser.parse_args()

o = Interbis(args.url, verify_certificates=False)
o.connect_to_datastore(args.login, args.password)

with open(settings_path, 'r') as file:
    default_settings = json.load(file)

settings_sample = o.get_sample("/ELN_SETTINGS/GENERAL_ELN_SETTINGS")
settings_sample.props["$eln_settings"] = json.dumps(default_settings)
settings_sample.save()

o.logout()
