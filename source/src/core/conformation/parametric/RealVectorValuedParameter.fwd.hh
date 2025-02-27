// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/parametric/RealVectorValuedParameter.fwd.hh
/// @brief  Owning pointers and whatnot for the class for holding a vector-valued parameter for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_conformation_parametric_RealVectorValuedParameter_fwd_hh
#define INCLUDED_core_conformation_parametric_RealVectorValuedParameter_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers

namespace core {
namespace conformation {
namespace parametric {

class RealVectorValuedParameter;

typedef  utility::pointer::weak_ptr< RealVectorValuedParameter >  RealVectorValuedParameterAP;
typedef  utility::pointer::weak_ptr< RealVectorValuedParameter const >  RealVectorValuedParameterCAP;
typedef  utility::pointer::shared_ptr< RealVectorValuedParameter >  RealVectorValuedParameterOP;
typedef  utility::pointer::shared_ptr< RealVectorValuedParameter const >  RealVectorValuedParameterCOP;

typedef  utility::vector1< RealVectorValuedParameterOP >  RealVectorValuedParameterOPs;
typedef  utility::vector1< RealVectorValuedParameterCOP >  RealVectorValuedParameterCOPs;
typedef  utility::vector1< RealVectorValuedParameterCAP >  RealVectorValuedParameterCAPs;

} // namespace parametric
} // namespace conformation
} // namespace core

#endif // INCLUDED_core_conformation_parametric_RealVectorValuedParameter_fwd_hh
