// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

# include <iostream>
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

int main( int argc, char ** argv ) {
    devel::init( argc, argv );
    utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
    if ( filenames.size() > 0 ) {
    std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
    } else {
        std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
        return 1;
    }

    // Read in pdb
    core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );

    // Set score funciton
    core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
    core::Real score = sfxn->score( *mypose );

    // Print score
    // mypose->energies().show( std::cout );
    std::cout << "Score Results is: " << score << std::endl;

    return 0;
}