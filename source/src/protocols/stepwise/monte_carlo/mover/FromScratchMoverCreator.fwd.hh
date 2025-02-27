// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/stepwise/monte_carlo/mover/FromScratchMoverCreator.fwd.hh
///
/// @brief
/// @author Andrew Watkins


#ifndef INCLUDED_protocols_stepwise_monte_carlo_mover_FromScratchMoverCreator_fwd_hh
#define INCLUDED_protocols_stepwise_monte_carlo_mover_FromScratchMoverCreator_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {


class FromScratchMoverCreator; // fwd declaration
typedef utility::pointer::shared_ptr< FromScratchMoverCreator > FromScratchMoverCreatorOP;
typedef utility::pointer::shared_ptr< FromScratchMoverCreator const > FromScratchMoverCreatorCOP;


} // namespace mover
} // namespace monte_carlo
} // namespace stepwise
} // namespace protocols

#endif // INCLUDED_protocols_stepwise_monte_carlo_mover_stepwise_monte_carlo_movers_FromScratchMover_FWD_HH
