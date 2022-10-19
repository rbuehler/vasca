"""
Tests validating the accessibility of configuration files
relevant to run the ResourceManager class.
"""
import os

import yaml

# paths
FILE_DIR = os.path.dirname(os.path.abspath(__file__))  # this file
TEST_DIR = FILE_DIR + "/.."  # test directory
PACKAGE_DIR = FILE_DIR + "/../.."  # root level of the repository


def test_global_paths():
    """
    Tests folder structure of the repository.
    """
    paths = [FILE_DIR, TEST_DIR, PACKAGE_DIR]
    assert all([os.path.isdir(dir) for dir in paths])


def test_load_env():
    """
    Tests if environment config file exists.
    """
    assert os.path.isfile(PACKAGE_DIR + "/.env")


def test_load_metadata():
    """
    Tests if the listed config files of the ResourceManager
    can be parsed with the YAML package and are successfully
    stored as sub-dictionaries of the metadata dictionary
    """
    # store metadata in dictionary
    # construct dict from list to keep it general
    # in case multiple metadata files are used in the future
    resources_names = ["catalog", "envs"]
    metadata = dict.fromkeys(resources_names)
    for key in metadata:
        with open(
            TEST_DIR + "/resource_metadata/resource_{}.yml".format(key), "r"
        ) as stream:
            try:
                metadata[key] = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print("YAMLerror:", exc)

    assert all([isinstance(metadata[key], dict) for key in metadata])
