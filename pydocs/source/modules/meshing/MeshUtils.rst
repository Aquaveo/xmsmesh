*********
MeshUtils
*********

The mesh_utils module contains the functions needed to generate a mesh. It also contains several
other helper functions that make using the meshing library easier. There are two ways to generate
a new mesh, :func:`mesh_utils.generate_2dm <xmsmesh.meshing.mesh_utils.generate_2dm>` and
:func:`mesh_utils.generate_mesh <xmsmesh.meshing.mesh_utils.generate_mesh>`.

The :func:`mesh_utils.generate_mesh <xmsmesh.meshing.mesh_utils.generate_mesh>` function takes a
:class:`MultiPolyMesherIo <xmsmesh.meshing.MultiPolyMesherIo>` object and generates the mesh by filling
in the :attr:`points <xmsmesh.meshing.MultiPolyMesherIo.points>` and the
:attr:`cells <xmsmesh.meshing.MultiPolyMesherIo.cells>` properties. The
:attr:`cell_polygons <xmsmesh.meshing.MultiPolyMesherIo.cell_polygons>` property is also generated if
you have set the :attr:`return_cell_polygons <xmsmesh.meshing.MultiPolyMesherIo.return_cell_polygons>`
property set to true on the MultiPolyMesherIo. This gives you
the definition of the mesh in memory so that you can work with it using other python libraries.

The :func:`mesh_utils.generate_mesh <xmsmesh.meshing.mesh_utils.generate_2dm>` function takes a
:class:`MultiPolyMesherIo <xmsmesh.meshing.MultiPolyMesherIo>` object and a file_name. While this
function does everything that generate_mesh does, it also generates a 2dm file which can be useful
if you are working with a product/library that can read 2dm files.

.. automodule:: xmsmesh.meshing.mesh_utils
   :members: