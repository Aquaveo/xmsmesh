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

// 4. External library headers

#include <boost/shared_ptr.hpp>

// 5. Shared code headers

#include <xmscore/misc/base_macros.h> // for XM_DISALLOW_COPY_AND_ASSIGN

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
class MeMultiPolyMesherIo;

//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------

//----- Function prototypes ----------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
class MeMultiPolyMesher
{
public:
  static boost::shared_ptr<MeMultiPolyMesher> New();
  /// \cond

  virtual bool MeshIt(MeMultiPolyMesherIo& a_io) = 0;
  virtual ~MeMultiPolyMesher() {}

  /// \endcond
protected:
  MeMultiPolyMesher() {}

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MeMultiPolyMesher)
}; // MeMultiPolyMesher

} // namespace xms
