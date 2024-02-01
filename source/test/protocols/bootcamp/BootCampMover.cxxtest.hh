// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/bootcamp/BootCampMover.cxxtest.hh
/// @brief  Mover to perform Monte Carlo relax based on DSSP SS fold-tree
/// @author Aseamann (austin.seamann@rutgers.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/kinematics/Edge.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>  // Needed for scoring
#include <core/scoring/ScoreFunction.hh> // Needed for scoring 

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/deep_copy.hh>

// C++ Headers
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/moves/MoverFactory.hh>

static basic::Tracer TR("BootCampMover");


class BootCampMover : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}

	void tearDown() {

	}

	// --------------- Test Cases --------------- //
	void test_bootcamp_mover() {
		// Setup bootcamp mover with MoverFactory
		protocols::moves::MoverOP base_mover_op = protocols::moves::MoverFactory::get_instance()->newMover("BootCampMover");

		protocols::bootcamp::BootCampMoverOP bcm_op = utility::pointer::dynamic_pointer_cast< protocols::bootcamp::BootCampMover > ( base_mover_op );

		// Confirm bcm_op is not NULL
		TS_ASSERT_DIFFERS( bcm_op, nullptr );
	}

	void test_bootcamp_mover_set_sfxn() {
		// Setup bootcamp mover with MoverFactory
		protocols::moves::MoverOP base_mover_op = protocols::moves::MoverFactory::get_instance()->newMover("BootCampMover");

		protocols::bootcamp::BootCampMoverOP bcm_op = utility::pointer::dynamic_pointer_cast< protocols::bootcamp::BootCampMover > ( base_mover_op );

		// Setup ScoreFunction
		core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();

		// Set ScoreFunction
		bcm_op->set_score_function( sfxn );

		// Confirm ScoreFunction is set
		TS_ASSERT_DIFFERS( bcm_op->get_score_function(), nullptr );
	}

	void test_bootcamp_mover_set_num_iter() {
		// Setup bootcamp mover with MoverFactory
		protocols::moves::MoverOP base_mover_op = protocols::moves::MoverFactory::get_instance()->newMover("BootCampMover");

		protocols::bootcamp::BootCampMoverOP bcm_op = utility::pointer::dynamic_pointer_cast< protocols::bootcamp::BootCampMover > ( base_mover_op );

		// Set num iterations
		bcm_op->set_num_iterations( 5 );

		// Confirm num iterations is set
		TS_ASSERT_EQUALS( bcm_op->get_num_iterations(), 5 );
	}



};
