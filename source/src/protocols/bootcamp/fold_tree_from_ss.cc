// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/bootcamp/fold_tree_from_ss.cc
/// @brief
/// @author Austin Seamann (austin.seamann@rutgers.edu)

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

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// Setup tracer cout
static basic::Tracer TR( "apps.pilot.aseamann" );

// Forward
namespace protocols {
namespace bootcamp {

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
    TR << "SS Element " << ii << " from residue "
      << ss_boundaries[ ii ].first << " to "
      << ss_boundaries[ ii ].second << std::endl;
  }
  return ss_boundaries;
}

core::kinematics::FoldTreeOP
fold_tree_from_dssp_string( std::string const & dssp_string )
{
	// Setup FoldTree
	core::kinematics::FoldTreeOP ft( new core::kinematics::FoldTree() );
	
	// Send string to identify_secondary_structure_spans
	utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries = identify_secondary_structure_spans( dssp_string );

	// Loop over pairs and add to fold tree
	utility::vector1< core::Size > w_loop_mid_points;
	for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
		// Get start and end of ss element
		core::Size start = ss_boundaries[ ii ].first;
		core::Size end = ss_boundaries[ ii ].second;
		core::Size mid = ( ( end - start ) / 2 ) + start;
		// Get the midpoints of the loops based on the last ss end and current ss start
		if ( ii > 1 ) {
			core::Size w_loop_mid = ( ( start - ss_boundaries[ ii - 1 ].second ) / 2 ) + ss_boundaries[ ii - 1 ].second;
			w_loop_mid_points.push_back( w_loop_mid );
		} 
		// Add mid point to fold tree
		w_loop_mid_points.push_back( mid );
	}

	// Build fold tree, starting at position 1 and to each mid point
	core::Size jump = 1;
	core::Size count_ss = 1;
	bool ss = true;
	core::Size current_mid;
	for ( core::Size ii = 1; ii <= w_loop_mid_points.size(); ++ii ) {

		current_mid = w_loop_mid_points[ ii ];

		// Check to see if current jump is within a ss or not, update bool ss
		if ( ii != 1 ) {
			if ( (ss_boundaries[count_ss - 1].first <= current_mid) && (current_mid >= ss_boundaries[count_ss - 1].second) ) {
				ss = true;
			} else {
				ss = false;
			}
		}

		if ( ii == 1 ) {
			// Add edge from 1 to first mid point
			ft->add_edge( current_mid, 1, core::kinematics::Edge::PEPTIDE );
			// Add edge from start node to end of ss element
			ft->add_edge( current_mid, ss_boundaries[ count_ss++ ].second, core::kinematics::Edge::PEPTIDE );
		} else if ( ss ) {
			// Add edge from first mid point to current mid point
			ft->add_edge( w_loop_mid_points[ 1 ], current_mid, jump++ );
			// Add current mid point to start/end of ss element
			ft->add_edge( current_mid, ss_boundaries[ count_ss - 1 ].second + 1 , core::kinematics::Edge::PEPTIDE );
			if (ii != w_loop_mid_points.size()) {
				ft->add_edge( current_mid, ss_boundaries[ count_ss ].first - 1 , core::kinematics::Edge::PEPTIDE );
			} else {
				ft->add_edge( current_mid, dssp_string.size(), core::kinematics::Edge::PEPTIDE );
			}
			// Update count_ss
			count_ss++;
		} else {
			// Add edge from first mid point to current mid point
			ft->add_edge( w_loop_mid_points[ 1 ], current_mid, jump++ );
			// Add current mid point to start/end of ss element
			ft->add_edge( current_mid, ss_boundaries[ count_ss - 1 ].first, core::kinematics::Edge::PEPTIDE );
			if (ii != w_loop_mid_points.size()) {
				ft->add_edge( current_mid, ss_boundaries[ count_ss - 1 ].second, core::kinematics::Edge::PEPTIDE );
			} else {
				ft->add_edge( current_mid, dssp_string.size(), core::kinematics::Edge::PEPTIDE );
			}
		}
	}

	return ft;
}

core::kinematics::FoldTreeOP
fold_tree_from_ss( core::pose::Pose const & pose )
{
	// Setup DSSP Call
	std::string dssp_str = core::scoring::dssp::Dssp( pose ).get_dssp_secstruct();
	// Call fold_tree_from_dssp_string
	core::kinematics::FoldTreeOP ft_dssp = fold_tree_from_dssp_string( dssp_str );

	return ft_dssp;
}

} //backrub
} //protocols