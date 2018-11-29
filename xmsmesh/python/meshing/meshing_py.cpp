//------------------------------------------------------------------------------
/// \file
/// \brief
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>
#include <xmscore/python/misc/PyUtils.h>
#include <xmsinterp/interpolate/InterpBase.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmsmesh/python/meshing/meshing_py.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------

void initMeshing(py::module &m) {
    initMeMeshUtils(m);
    initMeMultiPolyMesher(m);
    initMePolyInput(m);
    initMeRefinePoint(m);
    initMeMultiPolyMesherIo(m);
    initMeMultiPolyTo2dm(m);
    initMePolyMesher(m);
    initMePolyRedistributePts(m);
}
//------------------------------------------------------------------------------
/// \brief Create string from MeRefinePoint
/// \param[in] a_refinePoint: MeRefinePoint object
/// \return a string
//------------------------------------------------------------------------------
std::string PyReprStringFromMeRefinePoint(const xms::MeRefinePoint& a_refinePoint)
{
  std::string cmp_bool = a_refinePoint.m_createMeshPoint ? "True" : "False";
  std::stringstream ss;
  ss << "point: (" << a_refinePoint.m_pt.x << ", " << a_refinePoint.m_pt.y << ", " << a_refinePoint.m_pt.z << ")\n";
  ss << "size: " << a_refinePoint.m_size << "\n";
  ss << "create_mesh_point: " << cmp_bool;
  return ss.str();
} // PyReprStringFromMeRefinePoint
//------------------------------------------------------------------------------
/// \brief Create string from MePolyInput
/// \param[in] a_refinePoint: MePolyInput object
/// \return a string
//------------------------------------------------------------------------------
std::string PyReprStringFromMePolyInput(const xms::MePolyInput& a_polyInput)
{
  std::string szf = a_polyInput.m_sizeFunction == nullptr ? "none" : a_polyInput.m_sizeFunction->ToString();
  std::string elevf = a_polyInput.m_elevFunction == nullptr ? "none" : a_polyInput.m_elevFunction->ToString();
  std::stringstream ss;
  ss << "outside_poly: " << xms::StringFromVecPt3d(a_polyInput.m_outPoly);
  ss << "inside_poly: " << xms::StringFromVecPt3d2d(a_polyInput.m_insidePolys);
  ss << "bias: " << a_polyInput.m_bias << "\n";
  ss << "size_func: \n" << xms::PyReprStringFromInterpBase(a_polyInput.m_sizeFunction) << "\n";
  ss << "const_size_function: " << a_polyInput.m_constSizeFunction << "\n";
  ss << "const_size_bias: " << a_polyInput.m_constSizeBias << "\n";
  ss << "poly_corners: " << xms::StringFromVecInt(a_polyInput.m_polyCorners);
  ss << "elev_func: \n" << xms::PyReprStringFromInterpBase(a_polyInput.m_elevFunction) << "\n";
  ss << "bound_pts_to_remove: " << xms::StringFromVecPt3d(a_polyInput.m_boundPtsToRemove);
  ss << "remove_internal_four_triangle_pts: " << a_polyInput.m_removeInternalFourTrianglePts << "\n";
  ss << "poly_id: " << a_polyInput.m_polyId << "\n";
  ss << "seed_points: " << xms::StringFromVecPt3d(a_polyInput.m_seedPoints);
  ss << "relaxation_method: " << a_polyInput.m_relaxationMethod;
  return ss.str();
} // PyReprStringFromMePolyInput