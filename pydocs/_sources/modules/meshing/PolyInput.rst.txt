*********
PolyInput
*********

The PolyInput class contains information about polygons that you would like to mesh using the
:class:`MultiPolyMesherIo <xmsmesh.meshing.MultiPolyMesherIo>` class. This class defines the
polygon to be meshed. A polygon is made up of an outer polygon
(:attr:`outside_poly <xmsmesh.meshing.PolyInput.outside_poly>`), and a list of inner polygons
(:attr:`inside_polys <xmsmesh.meshing.PolyInput.inside_polys>`) representing
holes in the outer polygon. The other options available in this class help define
features of the polygons, such as :attr:`seed_poins <xmsmesh.meshing.PolyInput.seed_points>`,
:attr:`relaxation_method <xmsmesh.meshing.PolyInput.relaxation_method>`, and many more
that are described below.

.. autoclass:: xmsmesh.meshing.PolyInput
   :members:

   .. automethod:: __init__