#pragma once
//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 3. Standard library headers
#include <vector>

// 4. External library headers

#include <boost/shared_ptr.hpp>

// 5. Shared code headers

#include <xmscore/misc/base_macros.h>   // for XM_DISALLOW_COPY_AND_ASSIGN
#include <xmscore/misc/boost_defines.h> // BSHP

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
class Observer;
class MeMultiPolyMesherIo;

//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------

//----- Function prototypes ----------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// \brief Creates a 2dm file of a mesh from polygons
/// \see MeMultiPolyTo2dmImpl
class MeMultiPolyTo2dm
{
public:
  static boost::shared_ptr<MeMultiPolyTo2dm> New();

  /// \cond
  virtual bool Generate2dm(MeMultiPolyMesherIo& a_input, const std::string& a_outFileName,
                           int a_precision = 14) = 0;
  virtual bool Generate2dm(MeMultiPolyMesherIo& a_input, std::ostream& a_os,
                           int a_precision = 14) = 0;

  virtual void SetObserver(BSHP<Observer> a_) = 0;
  virtual ~MeMultiPolyTo2dm() {}

protected:
  MeMultiPolyTo2dm() {}

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MeMultiPolyTo2dm);
  /// \endcond
}; // MeMultiPolyTo2dm

} // namespace xms
