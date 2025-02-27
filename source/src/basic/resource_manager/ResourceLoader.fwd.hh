// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceLoader.fwd.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceLoader_FWD_HH
#define INCLUDED_basic_resource_manager_ResourceLoader_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace basic {
namespace resource_manager {

class ResourceLoader;
typedef utility::pointer::shared_ptr< ResourceLoader > ResourceLoaderOP;
typedef utility::pointer::shared_ptr< ResourceLoader const > ResourceLoaderCOP;


} // namespace resource_manager
} // namespace basic

#endif // INCLUDED_basic_resource_manager_ResourceLoader_FWD_HH
