#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from loguru import logger

from uvva import tables
from uvva.field import GALEXField
from uvva.tables import TableCollection


class Region(TableCollection):
    """
    :class: `~uvva.Region` defines a region in the sky as a
    list of uvva.field objects. It provides functionality to
    loop over fields to derive source lists, etc.
    """

    def __init__(self):
        """

        Notes
        -----
        Many class attributes are stored in astropy.table.Tables_. To see a
        description of each of their columns run :meth: `~uvva.Regions.info`.

        .. _astropy.table.Tables: https://docs.astropy.org/en/stable/api/astropy.table.Table.html

        Returns
        -------
        None.

        """
        # Sets skeleton
        super().__init__()

    @staticmethod
    def GALEXField_loader(field_id, obs):

        # Loads field according to load method specification
        # String matching is case insensitive (converts to all to lower case).
        if obs["field_options"]["load_method"].casefold() == "MAST".casefold():
            gf = GALEXField.from_MAST(
                obs_id=field_id,
                obs_filter=obs["obs_filter"],
                **obs["field_options"]["load_kwargs"],
            )
        elif obs["field_options"]["load_method"].casefold() == "UVVA".casefold():
            gf = GALEXField.from_UVVA(
                obs_id=field_id,
                obs_filter=obs["obs_filter"],
                **obs["field_options"]["load_kwargs"],
            )
        elif obs["field_options"]["load_method"].casefold() == "auto".casefold():
            pass
            # TODO: Lookahead via rm to check data availability.
            # Then "UVVA" is preferred for performance reasons.
            # Fallback to "MAST" & refresh=True if "UVVA" fails for some reason
            # (e.g. not complete set of tables stored in the fits file).
        else:
            raise ValueError(
                "Expected GALEXField load method specification from "
                "['MAST', 'UVVA', 'AUTO'], "
                f"got {obs['field_options']['load_method']}."
            )

        return gf

    @classmethod
    def load_from_config(cls, obs):
        """
        Loads region from configuration from dictionary

        Parameters
        ----------
        obs : dict
            Dictionary with region parameters derived from the uvva pipeline
            YAML configuration file.

        Returns
        -------
        None.

        """

        rg = cls()

        logger.debug("Loading fields from config file")

        if obs["observatory"] == "GALEX":
            rg.add_table(None, "region:tt_fields")
            rg.add_table(None, "region:tt_visits")

            # Temporary fix:
            obs["field_options"]["load_kwargs"]["load_products"] = False

            # Loop over fields and store info
            for field_id in obs["field_ids"]:

                gf = Region.GALEXField_loader(field_id, obs)
                field_info = dict(gf.tt_field[0])
                field_info["size"] = 0.55
                field_info["n_visits"] = gf.n_visits
                field_info["time_bin_size_sum"] = gf.time_bin_size_sum
                field_info["time_start"] = gf.time_start.mjd
                field_info["time_stop"] = gf.time_stop.mjd
                rg.tt_fields.add_row(field_info)

                # Loop over visits and store info
                keys_store = tables.dd_uvva_tables["base_field"]["tt_visits"]["names"]
                for ii in range(0, len(gf.tt_visits)):
                    visits_info = dict(gf.tt_visits[keys_store][0])
                    visits_info["field_id"] = field_id
                    rg.tt_visits.add_row(visits_info)
        else:
            logger.waring(
                "Selected observatory `" + obs["observatory"] + "` not supported"
            )

        return rg
