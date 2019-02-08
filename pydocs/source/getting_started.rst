Installation
------------

XmsMesh can be installed using `Anaconda <https://www.anaconda.com/download/>`_.

You can install XmsMesh using the `conda <https://www.anaconda.com/download/>`_ command::

   conda install -c aquaveo xmsmesh

This install XmsMesh and **all** the needed dependencies, including XmsCore, and NumPy.


Usage
-----

The meshing library takes MultiPolyMesherIos and converts them into meshes using the settings, and
the PolyInputs that are specified in the MultiPolyMesherIo.

The first step is to create a PolyInput that will contain your polygon definition. In this example
we will start with 4 corners and use the PolyRedistributePts object to fill in the spaces.

.. code-block:: python

    from xmsmesh.meshing import PolyInput
    from xmsmesh.meshing import PolyRedistributePts

    polygon_corners = [(0, 0, 0), (0, 100, 0), (100, 100, 0),
                       (0, 100, 0), (0, 0, 0)]
    rp = PolyRedistributePts()
    rp.set_constant_size_func(10)

    # We don't want a closed loop on our polygon_boundary
    # so we remove the last element
    polygon_boundary = rp.redistribute(polygon_corners)[:-1]

    poly_input = PolyInput(out_poly=polygon_boundary)

Now we will create our PolyMesherIo and generate our mesh

.. code-block:: python

    from xmsmesh.meshing import MultiPolyMesherIo
    from xmsmesh.meshing import mesh_utils

    mesher_io = MultiPolyMesherIo(poly_inputs=[poly_input])

    succeeded, errors = mesh_utils.generate_mesh(mesh_io=mesher_io)

In the above example **succeeded** is returned True if the meshing was successful, or False
if there were any errors. **errors** is a string errors that meshing returned. It is
an empty string if the meshing was successful.

The newly generated mesh can be accessed through the **points** property on the
MultiPolyMesher used to generate the mesh

.. code-block:: python

    print(mesher_io.points)
    print(mesher_io.cells)

There are many other options and settings that can be used to generate a mesh and the
documentation can be found in the **User Interface** section of this site. There are also
additional examples that can be found on the Examples_ page

.. _Examples: https://aquaveo.github.io/examples/xmsmesh/xmsmesh.html