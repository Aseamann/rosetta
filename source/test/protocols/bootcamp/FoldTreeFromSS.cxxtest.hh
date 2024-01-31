// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers
# include <utility/vector1.hh>

/// Project headers
#include <core/types.hh>

// C++ headers

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// -------------- Implemented Code --------------- //
utility::vector1< std::pair< core::Size, core::Size > >
identify_secondary_structure_spans( std::string const & ss_string )
{
  utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
  core::Size strand_start = -1;
  for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
    if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
      if ( int( strand_start ) == -1 ) {
        strand_start = ii;
      } else if ( ss_string[ii] != ss_string[strand_start] ) {
        ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
        strand_start = ii;
      }
    } else {
      if ( int( strand_start ) != -1 ) {
        ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
        strand_start = -1;
      }
    }
  }
  if ( int( strand_start ) != -1 ) {
    // last residue was part of a ss-eleemnt                                                                                                                                
    ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
  }
  for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
    std::cout << "SS Element " << ii << " from residue "
      << ss_boundaries[ ii ].first << " to "
      << ss_boundaries[ ii ].second << std::endl;
  }
  return ss_boundaries;
}


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

	}

};
