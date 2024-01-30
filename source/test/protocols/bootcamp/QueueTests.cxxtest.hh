// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  p/r/o/t/o/c/o/l/s///b/o/o/t/c/a/m/p/QueueTests.cxxtest.hh
/// @brief  Testing Queue for bootcamp
/// @author Aseamann (austin.seamann@rutgers.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/bootcamp/Queue.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <string>

static basic::Tracer TR("QueueTests");


class QueueTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();
		queue_ = protocols::bootcamp::Queue();
	}

	void tearDown() {

	}

	void test_first() {
		TS_TRACE( "Running my first unit test!" );
		TS_ASSERT( true );

	}

	void test_enqueue() {
		// Check if the first value adds
		std::string string1 = "Testing 1";
		queue_.enqueue(string1);
		TS_ASSERT_EQUALS(queue_.size(), 1);
		// Check if the second value adds
		std::string string2 = "Testing 2";
		queue_.enqueue(string2);
		TS_ASSERT_EQUALS(queue_.size(), 2);
	}

	void test_dequeue() {
		// Setup queue
		std::string string1 = "Testing 1";
		queue_.enqueue(string1);
		queue_.enqueue(string1);
		TS_ASSERT_EQUALS(queue_.size(), 2);
		// Check if we can dequeue the previously added strings
		queue_.dequeue();
		TS_ASSERT_EQUALS(queue_.size(), 1);
	}

	void test_is_empty() {
		// Setup queue
		std::string string1 = "Testing 1";
		queue_.enqueue(string1);
		// Check if is_empty is false
		TS_ASSERT(!queue_.is_empty());
		// Dequeue the single item
		queue_.dequeue();
		TS_ASSERT(queue_.is_empty());
	}

private:
	
	protocols::bootcamp::Queue queue_;

};
