import numpy as np
from matplotlib import colormaps as cm
from matplotlib.colors import hex2color
from matplotlib.lines import Line2D

import vasca.utils as vutils


def test_color_set():
    # Color palette for verification
    # fmt: off
    deep = [
        "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
        "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
    ]
    # fmt: on

    deep_rgba = [(*hex2color(c), 1.0) for c in deep]

    # Get set of colors
    cset_deep = vutils.color_palette("deep", len(deep))

    assert 4 == len(cset_deep[0])
    assert len(deep) == len(cset_deep)
    assert deep_rgba == cset_deep

    # Verify reversed color map
    cset_deep_r = vutils.color_palette("deep_r", len(deep))

    assert deep_rgba[-1] == cset_deep_r[0]

    # Verify with matplotlib viridis color codes
    c_viridis = cm.get_cmap("viridis").colors

    cset_viridis = vutils.color_palette("viridis", 256)

    assert tuple([*c_viridis[0], 1.0]) == cset_viridis[0]


def test_marker_set():
    # Test compatibility with matplotlib
    markers = vutils.marker_set(200)
    mpl_markers = Line2D.markers
    assert 200 == len(markers)
    assert all([m in mpl_markers for m in markers])

    # Test default exclusion
    assert all([m not in markers for m in [",", "8", "H"]])

    # Test custom exclusion
    exclude = ["*"]
    markers = vutils.marker_set(200, exclude=exclude, exclude_default=True)
    assert all([m not in markers for m in [",", "8", "H", "*"]])

    markers = vutils.marker_set(200, exclude=exclude, exclude_default=False)
    assert not all([m not in markers for m in [",", "8", "H", "*"]])


def test_utils_extr_value():
    # Test case 1: Finding the minimum value
    inputlist = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    assert vutils.extr_value(inputlist) == 1

    # Test case 2: Finding the maximum value
    inputlist = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    assert vutils.extr_value(inputlist, upper=True) == 9

    # Test case 3: Finding the minimum value with nested lists of different lengths
    inputlist = [[1, 2, 3], [4, 5], [6, 7, 8, 9]]
    assert vutils.extr_value(inputlist) == 1

    # Test case 4: Finding the maximum value with nested lists of different lengths
    inputlist = [[1, 2, 3], [4, 5], [6, 7, 8, 9]]
    assert vutils.extr_value(inputlist, upper=True) == 9

    # Test case 5: Empty input list should return None
    inputlist = []
    assert vutils.extr_value(inputlist) is None

    # Test case 6: Empty nested lists should return None
    inputlist = [[], []]
    assert vutils.extr_value(inputlist) is None

    # Test case 7: All nested lists contain NaN values, so the result should be None
    inputlist = [[np.nan, np.nan], [np.nan, np.nan]]
    assert vutils.extr_value(inputlist) is None

    # Test case 8: Custom input with a mix of positive, negative, and zero values
    inputlist = [[-1, 0, 1], [2, -3, 4], [5, -6, 7]]
    assert vutils.extr_value(inputlist) == -6

    # Test case 9: Custom input with a mix of positive, negative, and zero values,
    # finding the maximum
    inputlist = [[-1, 0, 1], [2, -3, 4], [5, -6, 7]]
    assert vutils.extr_value(inputlist, upper=True) == 7


def test_utils_get_hist_bins():
    # Test case 1: Testing with a list and default parameters
    data = [1.5, 2.5, 3.5, 4.5, 5.5]
    bin_size = 1
    bins = vutils.get_hist_bins(data, bin_size)
    assert np.array_equal(bins, np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))

    # Test case 2: Testing with an array and custom bin_scale
    data = np.array([0.25, 0.5, 0.75])
    bin_size = 0.1
    bins = vutils.get_hist_bins(data, bin_size)
    assert np.array_equal(bins, np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]))

    # Test case 3: Testing with a list, custom bin_size, and is_list=True
    data = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 12, 13]]
    bin_size = 2
    bins = vutils.get_hist_bins(data, bin_size)
    assert np.array_equal(bins, np.array([1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0]))


def test_utils_name2id():
    # Tests generic functionality (64-bit)
    name = "foo"
    id_int = vutils.name2id(name, bits=64)

    assert 64 >= id_int.bit_length()

    # Tests generic functionality (64-bit)
    name = "bar"
    id_int = vutils.name2id(name, bits=32)

    assert 32 >= id_int.bit_length()

    # Tests consistency (64-bit)
    name = "123"
    id_int = 11326344294570419622
    assert id_int == vutils.name2id(name, bits=64)

    # Tests consistency (32-bit)
    name = "123"
    id_int = 1503946150
    assert id_int == vutils.name2id(name, bits=32)

    # Tests error handling

    try:
        vutils.name2id("foo", bits=123)
    except ValueError as e:
        assert isinstance(e, ValueError)
