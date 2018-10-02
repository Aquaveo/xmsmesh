//------------------------------------------------------------------------------
/// \file
/// \brief
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <boost/shared_ptr.hpp>
#include <xmscore/python/misc/PyUtils.h>
#include <xmsgrid/ugrid/XmUGrid.h>
#include <xmsmesh/meshing/MeBadQuadRemover.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMeBadQuadRemover(py::module &m) {
  py::class_<xms::MeBadQuadRemover, boost::shared_ptr<xms::MeBadQuadRemover>>
    bad_quad_remover(m, "MeBadQuadRemover");

  bad_quad_remover.def(py::init<>([](boost::shared_ptr<xms::XmUGrid> ugrid)
  {
    return xms::MeBadQuadRemover::New(ugrid);
  }));

  // ---------------------------------------------------------------------------
  // function: remove_bad_quads
  // ---------------------------------------------------------------------------
  const char* remove_bad_quads_doc = R"pydoc(
    Remove bad quads and return a reconstructed UGrid with them removed.

    Args:
      max_aspect (double): The maximum aspect ratio for the diagonals to allowed
        for the quad to be collapsed.

    Returns:
      boost::shared_ptr<xms::XmUGrid> The reconstructed UGrid with the bad quads
        removed.

  )pydoc";
    bad_quad_remover.def("remove_bad_quads", [](xms::MeBadQuadRemover &self, double max_aspect = 0.7) -> boost::shared_ptr<xms::XmUGrid> {
      return self.RemoveBadQuads(max_aspect);
  }, remove_bad_quads_doc,py::arg("max_aspect"));
}