from pybis import Openbis
import os
from getpass import getpass

class LeBeOpenbis(Openbis):
    def __init__(self, url='https://test.datastore.bam.de/openbis/', verify_certificates=True, token=None, use_cache=True, allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks=False):
        super().__init__(url, verify_certificates, token, use_cache, allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks)
    
    def connect_to_datastore(self, *args, **kwargs):
        """Connects to a Datastore

        Args:
            url (str, optional): Datastore URL.
                    Defaults to 'https://test.datastore.bam.de/openbis/'.

        Returns:
            Openbis: Connected and running Datastore instance
        """

        # If a session is already active (when running a notebook or script again) return the active object
        # Instead of connecting again
        if not(self.is_session_active()):
            
            # We can parse the username and password as user input
            # If left empty the username will be grabbed from the os username
            # and password will have to be entered by the user
            if 'username' in kwargs:
                os.environ['OPENBIS_USERNAME'] = kwargs['username']
            else:
                os.environ['OPENBIS_USERNAME'] = os.getlogin()

            if 'password' in kwargs:
                os.environ['OPENBIS_PASSWORD'] = kwargs['password']
            else:
                os.environ['OPENBIS_PASSWORD'] = getpass("Give Password: ")

            try:
                self.login(os.environ['OPENBIS_USERNAME'],
                        os.environ['OPENBIS_PASSWORD'])
            except ValueError:
                print("Wrong Credentials")
                exit(1)

        return self
