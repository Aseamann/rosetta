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
# include <numeric/random/random.hh>  // For random number
# include <protocols/moves/MonteCarlo.hh>  // For montecarlo

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
    // Collect length of residues
    int N = mypose->size();

    // Set score funciton
    core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
    core::Real score = sfxn->score( *mypose );

    // Print score
    // mypose->energies().show( std::cout );
    std::cout << "Score Results is: " << score << std::endl;

    // Setup monte carlo object
    core::Real temp = 1.0;
    protocols::moves::MonteCarloOP mc ( new protocols::moves::MonteCarlo( *mypose, *sfxn, temp ) );

    int loop_count = 10;
    for( int i = 0 ; i < loop_count; i++ ) {
        // Determine random residue
        double uniform_random_number = numeric::random::uniform();
        core::Size randres = static_cast< core::Size > (uniform_random_number * N + 1 );
        // Collect random numbers to adjust phi/psi
        core::Real pert1 = numeric::random::gaussian();
        core::Real pert2 = numeric::random::gaussian();
        // Collect original phi/psi
        core::Real orig_phi = mypose->phi( randres );
        core::Real orig_psi = mypose->psi( randres );
        // Set random phi/psi
        mypose->set_phi( randres, orig_phi + pert1 );
        mypose->set_psi( randres, orig_psi + pert2 );

        mc->boltzmann( *mypose );
    }

    // Check pose
    core::Real score2 = sfxn->score( *mypose );
    std::cout << "Updated Score: " << score2 << std::endl;

    return 0;
}