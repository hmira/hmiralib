/*
 * iterable_vertices.h
 *
 *  Created on: Jul 24, 2012
 *      Author: hmirap
 */

#ifndef ITERABLE_VERTICES_H_
#define ITERABLE_VERTICES_H_

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

namespace hmira
{

namespace concepts
{

/*!
 *  \struct IterableVerticesConcept
 *  \brief concept required for a operation that iterates through all vertices
 * 
 *  implemented using boost::function_requires
 *  refines MeshConcept
 * 
 * \ingroup concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct IterableVerticesConcept
{
	typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
	typedef typename TMesh_Traits::vertex_iterator vertex_iterator;

	TMesh m;
	std::pair<vertex_iterator,vertex_iterator> vp;
	vertex_iterator vi;
	vertex_descriptor v;

	void constraints() {
		v = *vi;
		++vi;

		vp = TMesh_Traits::get_all_vertices(m);
	}
};

}//concepts

}//hmira

#endif /* ITERABLE_VERTICES_H_ */
