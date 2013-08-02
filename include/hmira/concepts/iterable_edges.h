/*
 * iterable_vertices.h
 *
 *  Created on: Jul 24, 2012
 *      Author: hmirap
 */

#ifndef ITERABLE_EDGES_H_
#define ITERABLE_EDGES_H_

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

namespace hmira
{

namespace concepts
{

/*!
 *  \struct IterableEdgesConcept
 *  \brief concept required for a operation that iterates through all edges
 * 
 *  implemented using boost::function_requires
 * 
 * \ingroup concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct IterableEdgesConcept
{
	typedef typename TMesh_Traits::edge_descriptor edge_descriptor;
	typedef typename TMesh_Traits::edge_iterator edge_iterator;

	TMesh m;
	std::pair<edge_iterator,edge_iterator> ep;
	edge_iterator ei;
	edge_descriptor e;

	void constraints() {
		e = *ei;
		++ei;
		ep = TMesh_Traits::get_all_edges(m);
	}
};

}//concepts

}//hmira

#endif /* ITERABLE_EDGES_H_ */
