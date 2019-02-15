*******************
PolyRedistributePts
*******************

The PolyRedistributePts class is used to help redistribute pts on a polyline.
There are several functions on this class used to set options
(:func:`set_constant_size_bias <xmsmesh.meshing.PolyRedistributePts.set_constant_size_bias>`,
:func:`set_constant_size_func <xmsmesh.meshing.PolyRedistributePts.set_constant_size_func>`,
:func:`set_size_func <xmsmesh.meshing.PolyRedistributePts.set_size_func>`) as well as the
function :func:`redistribute <xmsmesh.meshing.PolyRedistributePts.redistribute>` that actaully does all the work.

Example:

.. code-block:: python

    from xmsmesh.meshing import PolyRedistributePts

    polygon_corners = [(0, 0, 0), (0, 100, 0), (100, 100, 0),
                       (100, 0, 0), (0, 0, 0)]
    rp = PolyRedistributePts()
    rp.set_constant_size_func(10)
    redistributed_points = rp.redistribute(polygon_corners)


.. autoclass:: xmsmesh.meshing.PolyRedistributePts
   :members:
