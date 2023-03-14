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
