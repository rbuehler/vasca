"""
VASCA column definitions.

This module is currently not in productive use. Using this, the plan is to refactor the
tables dictionary module in the future.
"""

import dataclasses


@dataclasses.dataclass
class VASCAColumn:
    """Base class describing all parameters of VASCA columns.

    The parameters defined here, are used to construct the various :py:class:`~astropy.table.Table`
    objects to be held by a :class:`~vasca.tables.TableCollection`.

    :param name: Column Name
    :type name: str
    :param dtype: Data type of the column data (FITS compliant)
    :type dtype: str
    :param unit: Unit of the column data (FITS compliant)
    :type unit: str
    :param default: Default value taken as placeholder during initialization of
        :class:`~vasca.tables.TableCollection`.
    :type default: str
    :param description: Short description of the column's content
    :type description: str

    :example:
        >>> dummy_column = VASCAColumn(name="foo", dtype="int32", unit="1", default="-1", description="bar")
        >>> dummy_column.to_dict()
        {'name': 'foo',
         'dtype': 'int32',
         'unit': '1',
         'default': '-1',
         'description': 'bar'}
    """  # noqa: E501

    #: Column Name
    name: str
    #: Data Type (FITS compliant)
    dtype: str
    #: Unit (FITS compliant)
    unit: str
    #: Default value taken as placeholder during initialization of
    #: :class:`~vasca.tables.TableCollection`.
    default: str
    #: Short description of the column's content
    description: str

    def to_dict(self) -> dict[str, str]:
        """Returns a dictionary representation of the dataclass using
        :py:func:`dataclasses.asdict`.
        """
        return dataclasses.asdict(self)


field_id = VASCAColumn(
    name="field_id",
    dtype="S32",
    unit="1",
    default="none",
    description="Field source ID number",
)
"""
Field source ID number
"""
