#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Resource manager for VASCA
"""

import inspect
import os
from pprint import pprint

import dotenv
import yaml

# paths relative to the ResourceManager class
CLASS_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = CLASS_DIR + "/.."


class ResourceManager:
    """
    Manages access to external resources.

    Resources like catalog files or large metadata files may be stored on
    remote storage systems such as DESY Sync & Share (SAS) or LUSTRE.
    Users may synchronize these resources onto their local systems
    to perform software analysis on the data. To aid collaborative development,
    the `ResourceManager` provides an interface to find the locations of resources
    on the user's local system.

    Parameters
    ----------
    verbose : bool, optional
        Enable verbose printout to stdout (default is False for no printout).

    Attributes
    ----------
    flags : dict
        Dictionary holding all flags.
    metadata : dict
        Dictionary holding the parsed contents of the YAML-config files
        located at `CLASSDIR/resource_metadata`. Additional info my be added
        e.g. the set status of the env vars and the path they are pointing to.

    Notes
    -----
    The implementation requires a set of environment variables. The full list
    of variables is defined in `CLASSDIR/resource_metadata/resource_envs.yml`.
    Users need to make sure these variables are known to the system. This can be
    achieved in two ways. Either the user manually sets each variable via
    the shell config file (`.bashrc`, `zshrc`, etc.) or the instance of
    `ResourceManager` automatically sets the variables if the config file at
    `PACKAGEDIR/.env` is found.
    """

    def __init__(
        self,
        verbose=False,
    ):
        # parse arguments
        self.flags = dict()
        self.flags["verbose"] = verbose

        # load metadata, most importantly the resource catalog
        self.metadata = self._load_metadata()

        # set environment variables
        self._load_env()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # do some cleanup here
        pass

    def _load_metadata(self):
        # store metadata in dictionary
        # construct dict from list to keep it general
        # in case multiple metadata files are used in the future
        metadata = dict.fromkeys(["catalog", "envs"])
        for key in metadata:
            with open(
                CLASS_DIR + "/resource_metadata/resource_{}.yml".format(key), "r"
            ) as stream:
                try:
                    metadata[key] = yaml.safe_load(stream)
                except yaml.YAMLError as exc:
                    print("YAMLerror:", exc)

        if self.flags["verbose"]:
            pprint(metadata)

        return metadata

    def _load_env(self, config_path=PACKAGE_DIR + "/.env", overwrite=False):
        """
        Verify and set the required environment variables (env vars).

        ...
        """
        # log if env vars are already set
        self._log_env_status()
        env_status_list = [
            self.metadata["envs"][var]["set"] for var in self.metadata["envs"]
        ]
        # check if env var config file exists
        # If True: store listed env vars in dictionary
        found_config = True if os.path.isfile(config_path) else False
        if found_config:
            env_config = {**dotenv.dotenv_values(config_path)}

        # if all env vars were previously set, check if they match with env config
        # give warning if mismatch is detected
        # if no mismatch is occurs: nominal functionality of ResourceManager is assumed
        # and loading of env vars is not attempted
        if all(env_status_list) and found_config:
            is_mismatch = False
            for var in self.metadata["envs"]:
                if not os.environ.get(var) == env_config[var]:
                    is_mismatch = True
                    print(
                        str.format(
                            "WARNING(ResourceManager.{}): "
                            "Mismatch for environent variable '{}' "
                            "between the version that is set: \n{}\n"
                            "and the version listet in conf file: \n{}",
                            inspect.currentframe().f_code.co_name,
                            var,
                            os.environ.get(var),
                            env_config[var],
                        )
                    )
            if not is_mismatch:
                if self.flags["verbose"]:
                    print(
                        str.format(
                            "ResourceManager.{}: All env vars were previously set. "
                            "No mismatch with env var config file occured",
                            inspect.currentframe().f_code.co_name,
                        )
                    )
                return
        # if all env vars were previously set and no config file exists
        # nominal functionality of ResourceManager is assumed
        # and loading of env vars is not attempted
        elif all(env_status_list) and not found_config:
            if self.flags["verbose"]:
                print(
                    str.format(
                        "ResourceManager.{}: All env vars were previously set. "
                        "Missing env var config file.",
                        inspect.currentframe().f_code.co_name,
                    )
                )
            return
        # raise error if no env var was previously set and config file is also missing
        elif not all(env_status_list) and not found_config:
            raise Exception(
                str.format(
                    "Unable to set environemt variables {}. "
                    "Solve it by adding them to the config file at '{}' "
                    "or setting the variables manually in your shell config.",
                    [
                        var
                        for var in self.metadata["envs"]
                        if not self.metadata["envs"][var]["set"]
                    ],
                    config_path,
                )
            )
        # if not all env vars have been set previously but a config file exists
        # missing env vars are replaced
        elif not all(env_status_list) and found_config:
            if overwrite:
                dotenv.load_dotenv(config_path)
                self._log_env_status()
            else:
                missing_env_var_list = [
                    var
                    for var in self.metadata["envs"]
                    if not self.metadata["envs"][var]["set"]
                ]
                # raise error env var config is empty
                if len(list(env_config.keys())) == 0:
                    raise Exception(
                        str.format(
                            "Unable to set environemt variables {}. "
                            "Solve it by adding them to the config file at '{}' "
                            "or setting the variables manually in your shell config.",
                            [
                                var
                                for var in self.metadata["envs"]
                                if not self.metadata["envs"][var]["set"]
                            ],
                            config_path,
                        )
                    )
                for var in missing_env_var_list:
                    # raise error if missing env var is not found in config file
                    if var not in env_config:
                        raise Exception("Unable to set environemt variable '{}'.", var)
                    # set the missing env var
                    else:
                        os.environ[var] = env_config[var]
                self._log_env_status()

                # check for mismatch
                is_mismatch = False
                for var in self.metadata["envs"]:
                    if not os.environ.get(var) == env_config[var]:
                        is_mismatch = True
                        print(
                            str.format(
                                "WARNING(ResourceManager.{}): "
                                "Mismatch for environent variable '{}' "
                                "between the version that is set: \n{}\n"
                                "and the version listet in conf file: \n{}",
                                inspect.currentframe().f_code.co_name,
                                var,
                                os.environ.get(var),
                                env_config[var],
                            )
                        )

    def _log_env_status(self):
        """
        Verify and log if env vars are set.

        ...
        """
        for var in self.metadata["envs"]:
            is_set = False if os.environ.get(var) is None else True
            self.metadata["envs"][var]["set"] = is_set
            if self.flags["verbose"]:
                print("'{}' is set '{}'".format(var, is_set))
            self.metadata["envs"][var]["path"] = os.environ.get(var)

    def _check_resource_catalog(self):
        """
        Dummy function atm.

        Function to check the resource_catalog.yml file for consistency.
        Potential problems to look out for are:
            - duplicate resource IDs/names for a given storage
            - multiple resources with the same name at different storages
        Additionally check if cataloged files really exist
        """
        pass

    def _check_resoruce_env_info(self):
        """
        Verify if
            - not multiple env vars exist referring to one project and storage system
        """
        pass

    def get_path(self, resource, storage):
        """
        Returns the local path to a given resource at the corresponding storage system.

        The path is constructed from the local path to the synced remote directory
        and the resource path relative to that directory.

        Parameters
        ----------
        resource : str
            Name of the resource for which the path is requested.
        storage : str
            Name of the storage system where `resource` is saved at.

        Returns
        -------
        path : str
            Full local path to `resource`

        Raises
        ------
        KeyError
            When either one of resource or storage names is not
            listed in resource_catalog.yml
        ValueError
            When the inference of the needed environment variable fails
            or the resource could not be found at the inferred path.

        Notes
        -----
        This method requires specific environment variables to be declared
        on the users system during runtime.
        The variable names are listed in `CLASSDIR/resource_metadata/resource_envs.yml`.

        Examples
        --------
        >>> rm = ResourceManager()
        >>> rm.get_path(resource="gal_visits_list", storage="sas_cloud")
        str(/local/path/to/resource.file)

        """
        # validate storage name
        if storage not in self.metadata["catalog"]:
            raise KeyError(
                str.format(
                    ("Unknown storage system '{}'. Select one from {}."),
                    storage,
                    [strg for strg in list(self.metadata["catalog"].keys())],
                )
            )
        # list of all known resources: [<resource name>]
        resource_list = [
            self.metadata["catalog"][strg][id]["name"]
            for strg in self.metadata["catalog"]
            for id in self.metadata["catalog"][strg]
        ]
        # verbose list of resources: [<resource name>(<storage>:<resource ID>)]
        resource_list_verbose = [
            str.format(
                ("{}({}:{})"), self.metadata["catalog"][strg][id]["name"], strg, id
            )
            for strg in self.metadata["catalog"]
            for id in self.metadata["catalog"][strg]
        ]
        # validate resource name
        if resource not in resource_list:
            raise KeyError(
                str.format(
                    "Unknown resource '{}'. Select one from {}.",
                    resource,
                    resource_list_verbose,
                )
            )

        # get resource metadata
        resource_id = None
        resource_description = None
        resource_type = None
        project = None
        for id in self.metadata["catalog"][storage]:
            if resource == self.metadata["catalog"][storage][id]["name"]:
                project = self.metadata["catalog"][storage][id]["project"]
                resource_id = id
                resource_description = self.metadata["catalog"][storage][id][
                    "description"
                ]
                resource_type = self.metadata["catalog"][storage][id]["type"]
        if project is None:
            raise ValueError(
                str.format(
                    "No matching resource found for '{}' at storage system '{}'",
                    resource,
                    storage,
                )
            )
        # identify environment variable
        env_var = None
        for var in self.metadata["envs"]:
            if (storage == self.metadata["envs"][var]["storage"]) and (
                project == self.metadata["envs"][var]["project"]
            ):
                env_var = var
        if env_var is None:
            raise ValueError(
                str.format(
                    "No matching environment variable found "
                    "for resource '{}' of project '{}' at storage '{}'."
                    "Verify metadata files for consistency.",
                    resource,
                    project,
                    storage,
                )
            )
        # construct path
        path = (
            os.environ.get(env_var)
            + self.metadata["catalog"][storage][resource_id]["path"]
        )

        # validate path
        if resource_type == "file" and not os.path.isfile(path):
            raise ValueError("Not a file at '{}'".format(path))
        elif resource_type == "directory" and not os.path.isdir(path):
            raise ValueError("Not a directory at '{}'".format(path))

        # success message if verbose
        if self.flags["verbose"]:
            print(
                "Found a {} at '{}'. Resource description: {}".format(
                    resource_type, path, resource_description
                )
            )

        return path
