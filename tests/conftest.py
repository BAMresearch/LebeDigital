import pytest


def pytest_configure(config):
    # register an additional marker
    config.addinivalue_line(
        "markers", "login: mark test to only run when provided with login"
    )


"""
--login and --password command line arguments are username and password for
openBIS (that you use to log in normally)
"""


def pytest_addoption(parser):
    parser.addoption("--login", action="store", default="no_cl_login")
    parser.addoption("--password", action="store", default="no_cl_password")
    parser.addoption("--url", action="store", default="https://localhost:8443/openbis/")


# @pytest.fixture(scope='session', autouse=True)
# def pytest_generate_tests(metafunc):
#     # This is called for every test. Only get/set command line arguments
#     # if the argument is specified in the list of test "fixturenames".
#     option_value_login = metafunc.config.option.login
#     if 'login' in metafunc.fixturenames and option_value_login is not None:
#         metafunc.parametrize("login", [option_value_login])
#
#     option_value_password = metafunc.config.option.password
#     if 'password' in metafunc.fixturenames and option_value_password is not None:
#         metafunc.parametrize("password", [option_value_password])
