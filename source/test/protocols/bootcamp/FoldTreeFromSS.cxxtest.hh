// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/FoldTreeFromSS.cxxtest.hh
/// @brief
/// @author Austin Seamann (austin.seamann@rutgers.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>  // Output

/// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/kinematics/Edge.hh>

// C++ headers
#include <protocols/bootcamp/fold_tree_from_ss.hh>

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// Setup tracer cout
static basic::Tracer TR( "apps.pilot.aseamann" );

using namespace protocols::bootcamp;


// --------------- Test Class --------------- //

class FoldTreeFromSSTests : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	void test_hello_world() {
		TS_ASSERT( true );
	}

	void test_identify_secondary_structure_spans() {
		// Test 1
		// Setup object
		std::string ss_string = "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   ";
		utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries = identify_secondary_structure_spans( ss_string );

		// Validate num elements
		TS_ASSERT_EQUALS( ss_boundaries.size(), 7 );

		// Setup pairs and add to vector to loop over and compare
		utility::vector1< std::pair< core::Size, core::Size > > test_pairs;
		std::pair< core::Size, core::Size > pair1 = std::make_pair( 4, 8 );
		test_pairs.push_back( pair1 );
		std::pair< core::Size, core::Size > pair2 = std::make_pair( 12, 19 );
		test_pairs.push_back( pair2 );
		std::pair< core::Size, core::Size > pair3 = std::make_pair( 22, 26 );
		test_pairs.push_back( pair3 );
		std::pair< core::Size, core::Size > pair4 = std::make_pair( 36, 41 );
		test_pairs.push_back( pair4 );
		std::pair< core::Size, core::Size > pair5 = std::make_pair( 45, 55 );
		test_pairs.push_back( pair5 );
		std::pair< core::Size, core::Size > pair6 = std::make_pair( 58, 62 );
		test_pairs.push_back( pair6 );
		std::pair< core::Size, core::Size > pair7 = std::make_pair( 65, 68 );
		test_pairs.push_back( pair7 );

		// Loop over pairs and compare
		for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
			TS_ASSERT_EQUALS( ss_boundaries[ ii ].first, test_pairs[ ii ].first );
			TS_ASSERT_EQUALS( ss_boundaries[ ii ].second, test_pairs[ ii ].second );
		}

		// Test Two
		std::string ss_string2 = "HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH ";
		utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries2 = identify_secondary_structure_spans( ss_string2 );

		// Validate num elements
		TS_ASSERT_EQUALS( ss_boundaries2.size(), 7 );

		// Setup pairs and add to vector to loop over and compare
		utility::vector1< std::pair< core::Size, core::Size > > test_pairs2;
		std::pair< core::Size, core::Size > pair21 = std::make_pair( 1, 7 );
		test_pairs2.push_back( pair21 );
		std::pair< core::Size, core::Size > pair22 = std::make_pair( 11, 22 );
		test_pairs2.push_back( pair22 );
		std::pair< core::Size, core::Size > pair23 = std::make_pair( 29, 40 );
		test_pairs2.push_back( pair23 );
		std::pair< core::Size, core::Size > pair24 = std::make_pair( 41, 50 );
		test_pairs2.push_back( pair24 );
		std::pair< core::Size, core::Size > pair25 = std::make_pair( 51, 57 );
		test_pairs2.push_back( pair25 );
		std::pair< core::Size, core::Size > pair26 = std::make_pair( 59, 62 );
		test_pairs2.push_back( pair26 );
		std::pair< core::Size, core::Size > pair27 = std::make_pair( 63, 65 );
		test_pairs2.push_back( pair27 );

		// Loop over pairs and compare
		for ( core::Size ii = 1; ii <= ss_boundaries2.size(); ++ii ) {
			TS_ASSERT_EQUALS( ss_boundaries2[ ii ].first, test_pairs2[ ii ].first );
			TS_ASSERT_EQUALS( ss_boundaries2[ ii ].second, test_pairs2[ ii ].second );
		}

		// Test Three
		std::string ss_string3 = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE";
		utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries3 = identify_secondary_structure_spans( ss_string3 );

		// Validate num elements
		TS_ASSERT_EQUALS( ss_boundaries3.size(), 9 );

		// Setup pairs and add to vector to loop over and compare
		utility::vector1< std::pair< core::Size, core::Size > > test_pairs3;
		std::pair< core::Size, core::Size > pair31 = std::make_pair( 1, 9 );
		test_pairs3.push_back( pair31 );
		std::pair< core::Size, core::Size > pair32 = std::make_pair( 11, 18 );
		test_pairs3.push_back( pair32 );
		std::pair< core::Size, core::Size > pair33 = std::make_pair( 20, 28 );
		test_pairs3.push_back( pair33 );
		std::pair< core::Size, core::Size > pair34 = std::make_pair( 30, 30 );
		test_pairs3.push_back( pair34 );
		std::pair< core::Size, core::Size > pair35 = std::make_pair( 32, 36 );
		test_pairs3.push_back( pair35 );
		std::pair< core::Size, core::Size > pair36 = std::make_pair( 38, 38 );
		test_pairs3.push_back( pair36 );
		std::pair< core::Size, core::Size > pair37 = std::make_pair( 40, 40 );
		test_pairs3.push_back( pair37 );
		std::pair< core::Size, core::Size > pair38 = std::make_pair( 42, 42 );
		test_pairs3.push_back( pair38 );
		std::pair< core::Size, core::Size > pair39 = std::make_pair( 44, 51 );
		test_pairs3.push_back( pair39 );

		// Loop over pairs and compare
		for ( core::Size ii = 1; ii <= ss_boundaries3.size(); ++ii ) {
			TS_ASSERT_EQUALS( ss_boundaries3[ ii ].first, test_pairs3[ ii ].first );
			TS_ASSERT_EQUALS( ss_boundaries3[ ii ].second, test_pairs3[ ii ].second );
		}

	}

	// Test to perform running fold_tree_from_dssp_string
	void test_fold_tree_from_dssp_string() {
		// Setup string
		std::string ss_string = "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     ";
		// Setup fold tree
		core::kinematics::FoldTreeOP ft = fold_tree_from_dssp_string( ss_string );
		// Fold tree to string
		std::string ft_string = ft->to_string();
		TR << "Fold Tree Test Fold Tree From Dssp:\n" << ft_string << std::endl;
		// Validate fold tree, num edges
		TS_ASSERT_EQUALS( ft->size(), 38 );
	}

	// Test to perform running on test_in.pdb
	void test_fold_tree_from_dssp_on_test_in_pdb() {
		// Setup pose
		core::pose::Pose pose;
		pose = create_test_in_pdb_pose();

		// Setup string
		std::string ss_string = core::scoring::dssp::Dssp( pose ).get_dssp_secstruct();

		// Setup fold tree
		core::kinematics::FoldTreeOP ft = fold_tree_from_dssp_string( ss_string );

		// Check fold tree
		bool ft_check = ft->check_fold_tree();

	}

};
