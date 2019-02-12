*****************
MultiPolyMesherIo
*****************

The MultiPolyMesherIo contains does the options for meshing using XmsMesh library. This class uses
a collection of :class:`xmsmesh.meshing.PolyInput` objects, :class:`xmsmesh.meshing.RefinePoint` objects
and several other options that can be used to generate meshes.

This class is also used to retrieve a mesh after it has been generated. The points and cells
properties on this class are used to pass back the geometry of the mesh that is generated using
function from the :mod:`xmsmesh.meshing.mesh_utils` module.

The points property is just a list of (x, y, z) coordinates representing the points of the mesh. Cells,
is a cell stream that is used to define the mesh. The cell stream is a list of integers used to describe
the mesh.

For example a cell stream from a MultiPolyMesherIo might look something like this.


.. code-block:: bash

    >> print(mesh_io.cells)
    >> [5, 3, 0, 40, 1, 5, 3, 40, 0, 39, ...]

The cell stream is organized as follows: [<CELL_TYPE>, <NUMBER_OF_ELEMENTS> <ELEMENT__ID_1>, ... <ELEMENT_ID_N>,
<CELL_TYPE>, <NUMBER_OF_ELEMENTS> <ELEMENT_ID_1>, ... <ELEMENT_ID_N>, ...]. So our example above could be read
as CELL_TYPE = 5 (TRIANGLE), there are 3 elements, and their ID's are 0, 30, and 1, and so on.

The cell_polygons property is a list of ints referring to the index of the PolyInput that each cell was generated
from. The length of cell_polygons will be the total number of cells in your mesh. This property can be useful
because each mesh can be generated from many different polygons. It may sometimes be useful to know which polygon
each cell came from.

.. autoclass:: xmsmesh.meshing.MultiPolyMesherIo
   :members:

   .. automethod:: __init__