// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Basic includes
#include <iostream>
#include <numeric>

// Rosetta includes
# include <basic/Tracer.hh>  // Output
# include <basic/options/option.hh>
# include <utility/vector1.hh>
# include <basic/options/keys/in.OptionKeys.gen.hh>
# include <devel/init.hh>
# include <core/import_pose/import_pose.hh>  // importing pdbs
# include <utility/pointer/owning_ptr.hh>  // Needed for every time a smart pointer is used
# include <core/pose/Pose.hh>  // Needed for every time a smart pointer is used
# include <core/scoring/ScoreFunctionFactory.hh>  // Needed for scoring
# include <core/scoring/ScoreFunction.hh> // Needed for scoring 
# include <core/scoring/Energies.hh> // Showing energies from scoring
# include <numeric/random/random.hh>  // For random number
# include <protocols/moves/MonteCarlo.hh>  // For montecarlo
# include <protocols/moves/PyMOLMover.hh>  // PyMolObserver
# include <core/pack/task/PackerTask.hh>  // for setting up PackerTask
# include <core/pack/task/TaskFactory.hh>  // for setting up TaskFactory
# include <core/pack/pack_rotamers.hh>  // for packing
# include <core/kinematics/MoveMap.hh> // for minimizing
# include <core/optimization/AtomTreeMinimizer.hh> // for minimizing
# include <core/optimization/MinimizerOptions.hh> // for minimizing
# include <protocols/bootcamp/fold_tree_from_ss.hh>  // for fold_tree_from_ss
# include <core/pose/variant_util.hh>  // setup cutpoints
# include <core/scoring/ScoreTypes.hh>  // for setting weight
# include <protocols/jd2/JobDistributor.hh>  // for setting jd2

// Setup tracer cout
static basic::Tracer TR( "apps.pilot.aseamann" );

int main( int argc, char ** argv ) {
    devel::init( argc, argv );
    utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
    if ( filenames.size() > 0 ) {
    TR << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
    } else {
        TR << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
        return 1;
    }

    // Read in pdb
    core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );

    protocols::bootcamp::BootCampMoverOP bootcamp_mover( new protocols::bootcamp::BootCampMover );
    bootcamp_mover->apply( *mypose );

    protocols::jd2::JobDistributor::get_instance()->go(bootcamp_mover);

    return 0;
}