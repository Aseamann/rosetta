// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/MatchSet.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/MatchSet.hh>

// Project headers
#include <utility/graph/DisjointSets.hh>

// Utility headers
#include <utility/LexicographicalIterator.hh>
#include <utility/FixedSizeLexicographicalIterator.hh>
#include <utility/FixedSizeLexicographicalIterator.tmpl.hh>

#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>
#include <set>
#include <iostream>

namespace protocols {
namespace match {

/// @brief Non-member function that handles the logic for finding a nearby bin in 6D when
/// taking a step of a specified size from an existing bin.
core::Size
advance_to_neighbor_bin(
	Real6 orig_point,
	Real6 steps, // 0 if no step, non-zero if some real step
	numeric::geometry::hashing::SixDCoordinateBinner const & binner,
	numeric::geometry::hashing::Bin6D & next_bin
);


HitHasher::HitHasher() : initialized_( false ) {}
HitHasher::~HitHasher() = default;

void
HitHasher::set_bounding_box(
	BoundingBox const & bb
)
{
	debug_assert( ! initialized_ );
	bb_ = bb;
}

void
HitHasher::set_uniform_xyz_bin_width( Real bin_width )
{
	debug_assert( ! initialized_ );
	debug_assert( bin_width > 0 );
	for ( core::Size ii = 1; ii <= 3; ++ii ) { xyz_bin_widths_[ ii ] = bin_width; }
}

void
HitHasher::set_uniform_euler_angle_bin_width( Real bin_width_degrees )
{
	debug_assert( ! initialized_ );
	debug_assert( bin_width_degrees > 0 );
	for ( core::Size ii = 1; ii <= 3; ++ii ) { euler_bin_widths_[ ii ] = bin_width_degrees; }
}

void
HitHasher::set_xyz_bin_widths( Vector const & bin_widths )
{
	debug_assert( ! initialized_ );
	for ( core::Size ii = 1; ii <= 3; ++ii ) {
		debug_assert( bin_widths( ii ) > 0 );
		xyz_bin_widths_[ ii ] = bin_widths( ii );
	}

}

void
HitHasher::set_euler_bin_widths( Vector const & euler_bin_widths )
{
	debug_assert( ! initialized_ );
	for ( core::Size ii = 1; ii <= 3; ++ii ) {
		debug_assert( euler_bin_widths( ii ) > 0 );
		euler_bin_widths_[ ii ] = euler_bin_widths( ii );
	}
}

void
HitHasher::set_nhits_per_match( core::Size num_geometric_constraints )
{
	debug_assert( ! initialized_ );
	n_geometric_constraints_per_match_ = num_geometric_constraints;
}

void
HitHasher::initialize()
{
	debug_assert( ! initialized_ );

	initialized_ = true;

	hit_hashes_.resize( N_HASH_MAPS );

	/// Data for initializing the lower corner in xyz space
	Vector xyz_lower;
	Real3 half_xyzbin_widths = xyz_bin_widths_;
	for ( core::Size ii = 1; ii <= 3; ++ii ) { half_xyzbin_widths[ ii ] *= 0.5; }

	/// Data for initializing the lower corner in Euler space
	Size3 euler_offsets;
	utility::vector1< core::Size > six_twos( 6, 2 );
	utility::LexicographicalIterator lex( six_twos );

	Real6 bin_widths;
	bin_widths[ 1 ] = xyz_bin_widths_[ 1 ];   bin_widths[ 2 ] = xyz_bin_widths_[ 2 ];   bin_widths[ 3 ] = xyz_bin_widths_[ 3 ];
	bin_widths[ 4 ] = euler_bin_widths_[ 1 ]; bin_widths[ 5 ] = euler_bin_widths_[ 2 ]; bin_widths[ 6 ] = euler_bin_widths_[ 3 ];

	core::Size count = 1;

	//std::cout << "HitHasher initialize:";
	//for ( core::Size ii = 1; ii <= 6; ++ii ) std::cout << bin_widths[ ii ] << " "; std::cout << std::endl;

	while ( ! lex.at_end() ) {

		for ( core::Size ii = 1; ii <= 3; ++ii ) {
			xyz_lower( ii ) = bb_.lower()( ii ) + ( lex[ ii ] - 1 ) * ( half_xyzbin_widths[ ii ] );
		}
		for ( core::Size ii = 1, iip3 = 4; ii <= 3; ++ii, ++iip3 ) {
			euler_offsets[ ii ] = ( lex[ iip3 ] - 1 ); // 1 or 0 for offset or not.
		}

		//std::cout << "bounding box " << count << " lower: ";
		//for ( core::Size ii = 1; ii <= 3; ++ii ) std::cout << xyz_lower( ii ) << " ";
		//std::cout << "upper: ";
		//for ( core::Size ii = 1; ii <= 3; ++ii ) std::cout << bb_.upper()( ii ) << " ";
		//for ( core::Size ii = 1; ii <= 3; ++ii ) std::cout << euler_offsets[ ii ] << " "; std::cout << std::endl;

		BoundingBox bb( xyz_lower, bb_.upper() );

		hit_hashes_[ count ].first = utility::pointer::make_shared< numeric::geometry::hashing::SixDCoordinateBinner >( bb, euler_offsets, bin_widths );
		++lex;
		++count;
	}
}

void
HitHasher::insert_hit( core::Size geometric_constraint_id, Hit const * hit )
{
	runtime_assert( initialized_ );
	runtime_assert( hit );
	Real6 const & geom( hit->second() );
	Vector point( Vector( geom[ 1 ], geom[ 2 ], geom[ 3 ] ));
	if ( ! bb_.contains( point ) ) {
		utility_exit_with_message( "ERROR: Attempted to insert hit into HitHasher outside of the HitHasher bounding box!" );
	}

	for ( core::Size ii = 1; ii <= N_HASH_MAPS; ++ii ) {
		if ( hit_hashes_[ ii ].first->contains( geom ) ) {
			boost::uint64_t bin_index = hit_hashes_[ ii ].first->bin_index( geom );
			HitHash::iterator iter = hit_hashes_[ ii ].second.find( bin_index );
			if ( iter == hit_hashes_[ ii ].second.end() ) {

				/*if ( geometric_constraint_id != 1 ) {
				std::cout << "Weird -- this bin index " << bin_index << " should already have been inserted to hash " << ii;
				for ( core::Size jj = 1; jj <= 6; ++jj ) std::cout << " " << geom[ jj ]; std::cout << std::endl;
				} else {
				std::cout << "Inserting bin index: " << bin_index << " into hash " << ii << std::endl;
				}*/

				MatchSet ms( n_geometric_constraints_per_match_ );
				ms[ geometric_constraint_id ].push_back( hit );
				hit_hashes_[ ii ].second.insert( std::make_pair( bin_index, ms ));
			} else {
				iter->second[ geometric_constraint_id ].push_back( hit );
			}
		}
	}
}


/// @brief Insert a hits into a particular hash maps
void
HitHasher::insert_hit( core::Size which_hash_map, core::Size geometric_constraint_id, Hit const * hit )
{
	runtime_assert( initialized_ );
	runtime_assert( hit );
	Real6 const & geom( hit->second() );
	Vector point( Vector( geom[ 1 ], geom[ 2 ], geom[ 3 ] ));
	if ( ! bb_.contains( point ) ) {
		std::cerr << "Weird: Attempted to insert hit into HitHasher outside of the HitHasher bounding box!" << std::endl;
		std::cerr << "point: " << point.x() << " " << point.y() << " " << point.z() << std::endl;
		std::cerr << "bb: lower "<< bb_.lower().x() << " " << bb_.lower().y() << " " << bb_.lower().z() << std::endl;
		std::cerr << "bb: upper "<< bb_.upper().x() << " " << bb_.upper().y() << " " << bb_.upper().z() << std::endl;
		return;
		//utility_exit_with_message( "ERROR: Attempted to insert hit into HitHasher outside of the HitHasher bounding box!" );
	}

	if ( hit_hashes_[ which_hash_map ].first->contains( geom ) ) {
		boost::uint64_t bin_index = hit_hashes_[ which_hash_map ].first->bin_index( geom );
		HitHash::iterator iter = hit_hashes_[ which_hash_map ].second.find( bin_index );
		if ( iter == hit_hashes_[ which_hash_map ].second.end() ) {

			/*if ( geometric_constraint_id != 1 ) {
			std::cout << "Weird -- this bin index " << bin_index << " should already have been inserted to hash " << which_hash_map;
			for ( core::Size jj = 1; jj <= 6; ++jj ) std::cout << " " << geom[ jj ]; std::cout << std::endl;
			} else {
			std::cout << "Inserting bin index: " << bin_index << " into hash " << which_hash_map << std::endl;
			}*/

			MatchSet ms( n_geometric_constraints_per_match_ );
			ms[ geometric_constraint_id ].push_back( hit );
			hit_hashes_[ which_hash_map ].second.insert( std::make_pair( bin_index, ms ));
		} else {
			iter->second[ geometric_constraint_id ].push_back( hit );
		}
	}
}

void
HitHasher::clear_hash_map( core::Size which_hash_map )
{
	hit_hashes_[ which_hash_map ].second.clear();
}

HitHasher::HitHash::const_iterator
HitHasher::hit_hash_begin( core::Size which_hash_map ) const
{
	return hit_hashes_[ which_hash_map ].second.begin();
}

HitHasher::HitHash::const_iterator
HitHasher::hit_hash_end( core::Size which_hash_map ) const
{
	return hit_hashes_[ which_hash_map ].second.end();
}

//////////////////////////////////////////////////////////////////////
HitNeighborFinder::HitNeighborFinder() : initialized_( false ) {}
HitNeighborFinder::~HitNeighborFinder() = default;

void
HitNeighborFinder::set_bounding_box(
	BoundingBox const & bb
)
{
	debug_assert( ! initialized_ );
	bb_ = bb;
}

void
HitNeighborFinder::set_uniform_xyz_bin_width( Real bin_width )
{
	debug_assert( ! initialized_ );
	debug_assert( bin_width > 0 );
	for ( core::Size ii = 1; ii <= 3; ++ii ) { xyz_bin_widths_[ ii ] = bin_width; }
}

void
HitNeighborFinder::set_uniform_euler_angle_bin_width( Real bin_width_degrees )
{
	debug_assert( ! initialized_ );
	debug_assert( bin_width_degrees > 0 );
	for ( core::Size ii = 1; ii <= 3; ++ii ) { euler_bin_widths_[ ii ] = bin_width_degrees; }
}

void
HitNeighborFinder::set_xyz_bin_widths( Vector const & bin_widths )
{
	debug_assert( ! initialized_ );
	for ( core::Size ii = 1; ii <= 3; ++ii ) {
		debug_assert( bin_widths( ii ) > 0 );
		xyz_bin_widths_[ ii ] = bin_widths( ii );
	}

}

void
HitNeighborFinder::set_euler_bin_widths( Vector const & euler_bin_widths )
{
	debug_assert( ! initialized_ );
	for ( core::Size ii = 1; ii <= 3; ++ii ) {
		debug_assert( euler_bin_widths( ii ) > 0 );
		euler_bin_widths_[ ii ] = euler_bin_widths( ii );
	}
}

void HitNeighborFinder::add_hits( std::list< Hit > const & hitlist )
{
	debug_assert( all_hits_.empty() );
	debug_assert( initialized_ );

	core::Size count_hits = 0;
	for ( auto const & iter : hitlist ) {
		++count_hits;
		all_hits_.push_back( std::make_pair( count_hits, & iter ));
	}
	hash_hits();
}

void HitNeighborFinder::add_hits( HitPtrList const & /*hitptrlist*/ )
{}


void
HitNeighborFinder::initialize()
{
	debug_assert( ! initialized_ );

	initialized_ = true;

	/// Set all euler-offsets to 0 (no offset).
	Size3 euler_offsets( 0 );

	Real6 bin_widths;
	bin_widths[ 1 ] = xyz_bin_widths_[ 1 ];   bin_widths[ 2 ] = xyz_bin_widths_[ 2 ];   bin_widths[ 3 ] = xyz_bin_widths_[ 3 ];
	bin_widths[ 4 ] = euler_bin_widths_[ 1 ]; bin_widths[ 5 ] = euler_bin_widths_[ 2 ]; bin_widths[ 6 ] = euler_bin_widths_[ 3 ];

	binner_ = utility::pointer::make_shared< numeric::geometry::hashing::SixDCoordinateBinner >( bb_, euler_offsets, bin_widths );

}


utility::vector1< HitNeighborFinder::HitPtrList >
HitNeighborFinder::connected_components() const
{
	utility::graph::DisjointSets ds( all_hits_.size() );

	// we need to visit all 2^6=64 neighboring bins of the query halfbin; the lexicographical iterator
	// ensures we visit them all.
	Bin6D twos( 2 ); // numer of halfbin neighbors for each dimension
	utility::FixedSizeLexicographicalIterator< 6 > halfbin_lex( twos );
	Bin6D twopows( 1 );
	for ( core::Size ii = 5; ii >= 1; --ii ) twopows[ ii ] = twopows[ ii + 1 ] * 2;

	//typedef fixedsizearray1< boost::uint64_t, 2 > halfbin_index_touple; // pos1 = bin index; pos2 = halfbin subindex
	//typedef utility::OrderedTuple<

	utility::vector1< bool > visited_halfbins( 64, false );

	for ( HitHash::const_iterator iter = hash_.begin(), iter_end = hash_.end(); iter != iter_end; ++iter ) {
		HitIndexList const & hitlist = iter->second;
		if ( hitlist.empty() ) continue;
		std::fill( visited_halfbins.begin(), visited_halfbins.end(), false ); // mark all halfbins as not-yet visited.
		core::Size first_hit_index = hitlist.begin()->first;
		// Hit const * first_hit = hitlist.begin()->second; // Unused variable causes a warning.
		// Bin6D query_bin = binner_->bin6( first_hit->second() ); // Unused variable causes a warning.
		//boost::uint64_t bin_index = binner_->bin_index( query_bin );

		/// 1. Mark all the hits in this bin as members of the same CC.
		for ( auto const & hitlist_iter : hitlist ) {
			if ( ds.ds_find( first_hit_index ) == ds.ds_find( hitlist_iter.first ) ) continue;
			ds.ds_union( first_hit_index, hitlist_iter.first );
		}

		/// 2. Iterate across all hits, and for each halfbin, take one hit and look for all
		/// halfbin-neighbors of that hit.
		for ( auto const & hitlist_iter : hitlist ) {
			core::Size query_id = hitlist_iter.first;
			Hit const * query_hit = hitlist_iter.second;
			Bin6D query_halfbin = binner_->halfbin6( query_hit->second() );
			core::Size query_halfbin_index = 1;
			for ( core::Size ii = 1; ii <= 6; ++ii ) if ( query_halfbin[ ii ] == 1 ) query_halfbin_index += twopows[ ii ];

			// only explore the neighbor halfbins
			if ( visited_halfbins[ query_halfbin_index ] ) continue;
			visited_halfbins[ query_halfbin_index ] = true;

			Real6 halfsteps( binner_->halfbin_widths() );
			for ( core::Size ii = 1; ii <= 6; ++ii ) {
				if ( query_halfbin[ ii ] == 0 ) halfsteps[ ii ] *= -1;
			}

			halfbin_lex.begin();
			while ( !halfbin_lex.at_end() ) {
				Bin6D nbbin;
				core::Size skip_pos = find_next_bin( query_hit->second(), halfsteps, halfbin_lex, nbbin );
				if ( skip_pos != 0 ) {
					halfbin_lex.continue_at_dimension( skip_pos );
					continue;
				}

				/// now query the hash map and find all the hits in the neighbor bin (if any)
				boost::uint64_t nbindex = binner_->bin_index( nbbin );
				HitHash::const_iterator nbr_hits = hash_.find( nbindex );

				if ( nbr_hits != hash_.end() ) {
					// we have hits in the neighbor bin
					for ( auto const & nbr_hits_iter : nbr_hits->second ) {
						core::Size nbr_id = nbr_hits_iter.first;
						if ( ds.ds_find( query_id ) == ds.ds_find( nbr_id ) ) continue; // already neighbors
						Bin6D nb_halfbin = binner_->halfbin6( nbr_hits_iter.second->second() );
						if ( within_reach( query_halfbin, nb_halfbin, nbbin, halfbin_lex ) ) {
							ds.ds_union( query_id, nbr_id ); // mark these hits as neighbors
						}
					}
				}
				++halfbin_lex;

			} // end while ( !halfbin_lex.at_end() )


		}


	}

	/*for ( HitIndexList::const_iterator iter = all_hits_.begin(), iter_end = all_hits_.end();
	iter != iter_end; ++iter ) {
	core::Size query_id = iter->first;
	Hit const * query_hit = iter->second;
	Bin6D query_bin = binner_->bin6( query_hit->second() );
	Bin6D query_halfbin = binner_->halfbin6( query_hit->second() );
	Real6 halfsteps( binner_->halfbin_widths() );
	for ( core::Size ii = 1; ii <= 6; ++ii ) {
	if ( query_halfbin[ ii ] == 0 ) halfsteps[ ii ] *= -1;
	}

	halfbin_lex.begin();
	while ( !halfbin_lex.at_end() ) {
	Bin6D nbbin;
	core::Size skip_pos = find_next_bin( query_hit->second(), halfsteps, halfbin_lex, nbbin );
	if ( skip_pos != 0 ) {
	halfbin_lex.continue_at_dimension( skip_pos );
	continue;
	}

	/// now query the hash map and find all the hits in the neighbor bin (if any)
	boost::uint64_t nbindex = binner_->bin_index( nbbin );
	HitHash::const_iterator nbr_hits = hash_.find( nbindex );

	if ( nbr_hits != hash_.end() ) {
	// we have hits in the neighbor bin
	for ( HitIndexList::const_iterator nbr_hits_iter = nbr_hits->second.begin(),
	nbr_hits_iter_end = nbr_hits->second.end();
	nbr_hits_iter != nbr_hits_iter_end; ++nbr_hits_iter ) {
	core::Size nbr_id = nbr_hits_iter->first;
	if ( ds.ds_find( query_id ) == ds.ds_find( nbr_id ) ) continue; // already neighbors
	Bin6D nb_halfbin = binner_->halfbin6( nbr_hits_iter->second->second() );
	if ( within_reach( query_halfbin, nb_halfbin, nbbin, halfbin_lex )) {
	ds.ds_union( query_id, nbr_id ); // mark these hits as neighbors
	}
	}
	}
	++halfbin_lex;

	} // end while ( !halfbin_lex.at_end() )
	}*/

	/// Iterate across all bins.  Then iterate across all halfbins, skipping over halfbins that have already been visited.

	// OK: now how many sets do we have, and who are their representatives?
	core::Size const nccs = ds.n_disjoint_sets(); // number of connected-components.
	utility::vector1< core::Size > ds_representatives_to_setnos( ds.n_nodes(), 0 );
	core::Size count_set = 0;
	for ( core::Size ii = 1; ii <= ds.n_nodes(); ++ii ) {
		if ( ds.ds_find( ii ) == ii ) {
			++count_set;
			ds_representatives_to_setnos[ ii ] = count_set;
		}
	}

	utility::vector1< HitPtrList > hit_ccs( nccs );
	for ( auto const & all_hit : all_hits_ ) {
		core::Size iterrep = ds.ds_find( all_hit.first );
		core::Size setid = ds_representatives_to_setnos[ iterrep ];
		debug_assert( setid != 0 );
		hit_ccs[ setid ].push_back( all_hit.second );
	}
	return hit_ccs;
}


/// @brief Find the neighbors of the given set of query hits.  This search iterates
/// across both the upper and the lower neighbors of the query hits (3^6 neighbors).
HitNeighborFinder::HitPtrList
HitNeighborFinder::neighbor_hits( HitPtrList const & queryhits ) const
{
	HitPtrList neighbors;
	std::set< core::Size > neighbor_set;
	std::set< boost::uint64_t > exhausted_bins;
	// we need to visit all 2^6=64 neighboring bins of the query halfbin; the lexicographical iterator
	// ensures we visit them all.
	Bin6D twos( 2 ); // numer of halfbin neighbors for each dimension
	utility::FixedSizeLexicographicalIterator< 6 > halfbin_lex( twos );

	for ( auto query_hit : queryhits ) {
		//core::Size query_id = iter->first;
		// Bin6D query_bin = binner_->bin6( query_hit->second() ); // Unused variable causes a warning.
		Bin6D query_halfbin = binner_->halfbin6( query_hit->second() );
		Real6 halfsteps( binner_->halfbin_widths() );
		for ( core::Size ii = 1; ii <= 6; ++ii ) {
			if ( query_halfbin[ ii ] == 0 ) halfsteps[ ii ] *= -1;
		}

		halfbin_lex.begin();
		while ( !halfbin_lex.at_end() ) {
			Bin6D nbbin;
			core::Size skip_pos = find_next_bin( query_hit->second(), halfsteps, halfbin_lex, nbbin );
			if ( skip_pos != 0 ) {
				halfbin_lex.continue_at_dimension( skip_pos );
				continue;
			}

			/// now query the hash map and find all the hits in the neighbor bin (if any)
			boost::uint64_t nbindex = binner_->bin_index( nbbin );
			HitHash::const_iterator nbr_hits = hash_.find( nbindex );

			if ( exhausted_bins.find( nbindex ) != exhausted_bins.end() ) { ++halfbin_lex; continue; }

			if ( nbr_hits != hash_.end() ) {
				// we have hits in the neighbor bin
				bool all_inserted = true;
				for ( auto const & nbr_hits_iter : nbr_hits->second ) {
					core::Size nbr_id = nbr_hits_iter.first;
					if ( neighbor_set.find( nbr_id ) != neighbor_set.end() ) {
						// already found this neighbor; note does not invalidate the "all_inserted" boolean
						// since this hit need not be examined again in the future.
						continue;
					}

					Bin6D nb_halfbin = binner_->halfbin6( nbr_hits_iter.second->second() );
					if ( within_reach( query_halfbin, nb_halfbin, nbbin, halfbin_lex ) ) {
						neighbor_set.insert( nbr_id );
						neighbors.push_back( nbr_hits_iter.second );
					} else {
						all_inserted = false; // a hit is out of range; don't add it.
					}
				}
				if ( all_inserted ) {
					// mark this bin as having no hits that need to be further examined
					exhausted_bins.insert( nbindex );
				}
			} else {
				/// If this bin is empty, we may as well avoid the hash lookup.
				exhausted_bins.insert( nbindex );
			}
			++halfbin_lex;

		} // end while ( !halfbin_lex.at_end() )
	}

	return neighbors;
}

void
HitNeighborFinder::hash_hits()
{
	for ( HitIndexList::const_iterator iter = all_hits_.begin(), iter_end = all_hits_.end();
			iter != iter_end; ++iter ) {
		boost::uint64_t bin_index = binner_->bin_index( iter->second->second() );
		HitHash::iterator hash_iter = hash_.find( bin_index );
		if ( hash_iter == hash_.end() ) {
			HitIndexList hitlist;
			hitlist.push_back( *iter );
			hash_[ bin_index ] = hitlist;
		} else {
			hash_iter->second.push_back( *iter );
		}
	}

}

/// @brief Find the neighbor bin for a point by taking a step along a certain offset vector
/// for a subset of the dimensions as given by a lex iterator (the halfbin_lex).
/// Returns 0 if successful, and the index of the euclidean dimension which is out-of-bounds
/// if the step is not successful.
HitNeighborFinder::Size
HitNeighborFinder::find_next_bin(
	Real6 orig_point,
	Real6 offsets,
	utility::FixedSizeLexicographicalIterator< 6 > const & halfbin_lex,
	Bin6D & next_bin
) const
{
	Real6 steps( 0.0 );
	for ( core::Size ii = 1; ii <= 6; ++ii ) {
		if ( halfbin_lex[ ii ] == 2 ) steps[ ii ] = offsets[ ii ];
	}
	return advance_to_neighbor_bin( orig_point, steps, *binner_, next_bin );

}


bool
HitNeighborFinder::within_reach(
	Bin6D const & query_halfbin,
	Bin6D const & nb_halfbin,
	Bin6D const & nbbin,
	utility::FixedSizeLexicographicalIterator< 6 > const & halfbin_lex
) const
{
	core::Size const NEIGHBOR_BIN = 2; // 1 will denote "stay in the same bin; 2 will denote look in the neighbor bin"

	// check everything but theta!
	for ( core::Size ii = 1; ii <= 5; ++ii ) {
		if ( halfbin_lex[ ii ] == NEIGHBOR_BIN && query_halfbin[ ii ] == nb_halfbin[ ii ] ) {
			/// Then this neighbor hit is too far away from the query hit: we've spanned
			/// the bin boundary to a neighboring bin, but the neighbor hit is on
			/// the other side of the halfbin divide from the query hit -- either
			/// halfbin[ ii ] == 0 && nbhalfbin[ ii ]
			return false;
		}
	}

	// Theta comparison: worry about wrapping!
	if ( halfbin_lex[ 6 ] == NEIGHBOR_BIN ) {
		if ( nbbin[ 6 ] == 0 || nbbin[ 6 ] + 1 == binner_->dimsizes()[ 6 ] ) {
			// theta has wrapped, so we're only within range as long as both halfbins agree
			if ( query_halfbin[ 6 ] != nb_halfbin[ 6 ] ) {
				return false;
			}
		} else {
			// theta did not wrap, so the two halfbins need to disagree to be within range.
			if ( query_halfbin[ 6 ] == nb_halfbin[ 6 ] ) {
				return false;
			}
		}
	}
	return true;
}


//////////////////////////////////////////////////////////////////////
MatchCounter::MatchCounter() : n_geom_csts_( 0 ), initialized_( false ) {}
MatchCounter::~MatchCounter() = default;

void
MatchCounter::set_bounding_box(
	BoundingBox const & bb
)
{
	debug_assert( ! initialized_ );
	bb_ = bb;
}

void
MatchCounter::set_n_geometric_constraints( core::Size ngeomcsts )
{
	debug_assert( ! initialized_ );
	debug_assert( n_geom_csts_ == 0 ); // this function should only be called once
	n_geom_csts_ = ngeomcsts;
}

void
MatchCounter::set_uniform_xyz_bin_width( Real bin_width )
{
	debug_assert( ! initialized_ );
	debug_assert( bin_width > 0 );
	for ( core::Size ii = 1; ii <= 3; ++ii ) { xyz_bin_widths_[ ii ] = bin_width; }
}

void
MatchCounter::set_uniform_euler_angle_bin_width( Real bin_width_degrees )
{
	debug_assert( ! initialized_ );
	debug_assert( bin_width_degrees > 0 );
	for ( core::Size ii = 1; ii <= 3; ++ii ) { euler_bin_widths_[ ii ] = bin_width_degrees; }
}

void
MatchCounter::set_xyz_bin_widths( Vector const & bin_widths )
{
	debug_assert( ! initialized_ );
	for ( core::Size ii = 1; ii <= 3; ++ii ) {
		debug_assert( bin_widths( ii ) > 0 );
		xyz_bin_widths_[ ii ] = bin_widths( ii );
	}

}

void
MatchCounter::set_euler_bin_widths( Vector const & euler_bin_widths )
{
	debug_assert( ! initialized_ );
	for ( core::Size ii = 1; ii <= 3; ++ii ) {
		debug_assert( euler_bin_widths( ii ) > 0 );
		euler_bin_widths_[ ii ] = euler_bin_widths( ii );
	}
}

void MatchCounter::add_hits( core::Size geomcst_id, std::list< Hit > const & hitlist )
{
	debug_assert( initialized_ );

	for ( auto const & iter : hitlist ) {
		boost::uint64_t bin_index = binner_->bin_index( iter.second() );
		HitHash::iterator hash_iter = hash_.find( bin_index );
		if ( hash_iter == hash_.end() ) {
			HitCounts & hitcount_v = hash_[ bin_index ];
			hitcount_v.resize( n_geom_csts_ );
			std::fill( hitcount_v.begin(), hitcount_v.end(), 0 );
			hitcount_v[ geomcst_id ] = 1;
		} else {
			hash_iter->second[ geomcst_id ] += 1;
		}
	}
}

void MatchCounter::add_hits( core::Size geomcst_id, std::list< Hit const * > const & hitlist )
{
	debug_assert( initialized_ );

	for ( auto iter : hitlist ) {
		boost::uint64_t bin_index = binner_->bin_index( iter->second() );
		HitHash::iterator hash_iter = hash_.find( bin_index );
		if ( hash_iter == hash_.end() ) {
			HitCounts & hitcount_v = hash_[ bin_index ];
			hitcount_v.resize( n_geom_csts_ );
			std::fill( hitcount_v.begin(), hitcount_v.end(), 0 );
			hitcount_v[ geomcst_id ] = 1;
		} else {
			hash_iter->second[ geomcst_id ] += 1;
		}
	}
}

/// @details Construct the hash table based on the parameters set previously
void
MatchCounter::initialize()
{
	debug_assert( ! initialized_ );

	initialized_ = true;

	/// Set all euler-offsets to 0 (no offset).
	Size3 euler_offsets( 0 );

	Real6 bin_widths;
	bin_widths[ 1 ] = xyz_bin_widths_[ 1 ];   bin_widths[ 2 ] = xyz_bin_widths_[ 2 ];   bin_widths[ 3 ] = xyz_bin_widths_[ 3 ];
	bin_widths[ 4 ] = euler_bin_widths_[ 1 ]; bin_widths[ 5 ] = euler_bin_widths_[ 2 ]; bin_widths[ 6 ] = euler_bin_widths_[ 3 ];

	binner_ = utility::pointer::make_shared< numeric::geometry::hashing::SixDCoordinateBinner >( bb_, euler_offsets, bin_widths );

}

/// @details This function computes a lower bound on the number of matches. It will miss matches
/// that are neighboring grid cells -- matches that the match enumeration step later will find.
/// This function, though, needs to be fast more than it needs to be perfectly accurate. A
/// previous implementation performed a perfect, non-redundant counting of all the matches,
/// but at considerable computational expense: at each decrease in the grid spacing, the
/// number of non-empty cells in the hash table increases. The result was that the dynamic-grid
/// refinement step would slow more and more as finer grid spacings were considered: exactly
/// opposite of the purpose of the dynamic grid refinment.
MatchCounter::Size
MatchCounter::count_n_matches() const
{
	debug_assert( initialized_ );
	core::Size const two_billion = 2000000000;
	Real const ln2e9 = 21.4; // e^21.4 ~ 2 billion

	core::Size total = 0;
	core::Size last_total = 0;
	for ( HitHash::const_iterator iter = hash_.begin(), iter_end = hash_.end(); iter != iter_end; ++iter ) {
		HitCounts const & hit_counts = iter->second;
		core::Size counts_this_bin = 1;
		Real log_counts_this_bin = 0;
		for ( core::Size ii = 1; ii <= hit_counts.size(); ++ii ) {
			counts_this_bin *= hit_counts[ii];
			if ( hit_counts[ii] != 0 ) {
				log_counts_this_bin += log((Real)hit_counts[ii]);
			}
		}
		if ( log_counts_this_bin > ln2e9 ) {
			return two_billion;
		}
		total += counts_this_bin;
		if ( total < last_total ) {
			// then we have had an integer overflow
			return two_billion;
		}
		last_total = total;
	}
	return total;

}

////////////////


core::Size
advance_to_neighbor_bin(
	Real6 orig_point,
	Real6 steps, // 0 if no step, non-zero if some real step
	numeric::geometry::hashing::SixDCoordinateBinner const & binner,
	numeric::geometry::hashing::Bin6D & next_bin
)
{
	using namespace core;

	numeric::geometry::hashing::Real6 alt_point( orig_point );
	numeric::geometry::hashing::Real3 orig_euler( 0.0 );
	numeric::geometry::hashing::Real3 euler_offsets( 0.0 );
	for ( core::Size ii = 1; ii <= 3; ++ii ) {
		alt_point[ ii ] += steps[ ii ];
		euler_offsets[ ii ] = steps[ ii + 3 ];
		orig_euler[ ii ] = orig_point[ ii + 3 ];
	}

	/// are we in range? -- return out-of-range index for advancing the lex if so
	for ( core::Size ii = 1; ii <= 3; ++ii ) {
		if ( binner.bounding_box().lower()( ii ) > alt_point[ ii ] ||
				binner.bounding_box().upper()( ii ) < alt_point[ ii ] ) {
			return ii;
		}
	}
	/// ok -- we're good!

	/// Complicated logic in advancing the euler angles; defer to another function
	numeric::geometry::hashing::Real3 new_euler = advance_euler_angles( orig_euler, euler_offsets );
	for ( core::Size ii = 1; ii <= 3; ++ii ) {
		alt_point[ ii + 3 ] = new_euler[ ii ];
	}

	next_bin = binner.bin6( alt_point );

	return 0;

}

/// @details "Advance" the euler angles by a given offset (postive or negative), and wrap
/// them back into their appropriate ranges ([0,360) for phi/psi, [0..180] for theta)
/// as necessary.
/// Requirements: orig_angles[ i ] + offsets[ i ] < 520 (360) && > -360 (-180) for phi/psi (theta).
numeric::geometry::hashing::Real3
advance_euler_angles(
	numeric::geometry::hashing::Real3 const & orig_angles,
	numeric::geometry::hashing::Real3 const & offsets
)
{
	using core::Size;

	numeric::geometry::hashing::Real3 new_euler_angles( orig_angles );
	for ( core::Size ii = 1; ii <= 3; ++ii ) new_euler_angles[ ii ] += offsets[ ii ];

	if ( new_euler_angles[ 3 ] < 0 || new_euler_angles[ 3 ] > 180 ) {
		/// 1st handle theta wrapping.
		if ( new_euler_angles[ 3 ] < 0 ) {
			new_euler_angles[ 3 ] *= -1.0; // wrap the angle back to positive values
		} else { // new_euler_angles[ 3 ] > 180
			debug_assert( new_euler_angles[ 3 ] < 360 );
			new_euler_angles[ 3 ] = 360 - new_euler_angles[ 3 ]; // 182 wraps to 178...
		}
		/// wrap phi/psi
		for ( core::Size ii = 1; ii <=2; ++ii ) {
			if ( new_euler_angles[ ii ] < 180 ) {
				new_euler_angles[ ii ] += 180;
			} else {
				new_euler_angles[ ii ] -= 180;
			}
		}
	} else {
		/// Theta doesn't wrap, so make sure phi and psi are in their proper ranges.
		for ( core::Size ii = 1; ii <= 2; ++ii ) {
			if ( new_euler_angles[ ii ] > 360 ) {
				new_euler_angles[ ii ] -= 360;
			} else if ( new_euler_angles[ ii ] < 0 ) {
				new_euler_angles[ ii ] += 360;
			}
		}
	}

	//std::cout << "orig_angles:"; for ( core::Size ii = 1; ii <= 3; ++ii ) std::cout << " " << orig_angles[ ii ]; std::cout << std::endl;
	//std::cout << "offsets:"; for ( core::Size ii = 1; ii <= 3; ++ii ) std::cout << " " << offsets[ ii ]; std::cout << std::endl;
	//std::cout << "new_euler_angles:"; for ( core::Size ii = 1; ii <= 3; ++ii ) std::cout << " " << new_euler_angles[ ii ]; std::cout << std::endl;
	return new_euler_angles;
}

////////////////////////////////////////////////////////////////////////////////////


MatchOutputTracker::MatchOutputTracker() = default;

void
MatchOutputTracker::note_output_match( match_lite const & m )
{
	MatchHash::const_iterator iter = hash_.find( m );
	if ( iter == hash_.end() ) {
		hash_.insert( std::make_pair( m, true ) );
	}
}

bool
MatchOutputTracker::match_has_been_output( match_lite const & m ) const
{
	return hash_.find( m ) != hash_.end();
}


}
}
