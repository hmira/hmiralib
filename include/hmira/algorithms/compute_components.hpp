/*
 * compute_components.h
 *
 *  Created on: Jul 20, 2012
 *      Author: hmirap
 */

#ifndef __COMPUTE_COMPONENTS_HPP__
#define __COMPUTE_COMPONENTS_HPP__

#include <vector>
#include <deque>
#include <algorithm>

#include <hmira/abstract/mesh_traits.hpp>
#include <hmira/concepts/algorithms/compute_components.hpp>

/*
 * Simple algorithm that computes the count of components of
 * the Mesh. Isolated vertex is considered as a component
 * Algorithm works as modified BFS
 */

namespace hmira
{

namespace algorithms
{

/**
 * 
 * \ingroup algorithms
 * 
 * template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>> int compute_components(TMesh& m)
 * \brief simple algorithm that computes the count of components
 * \param m a polygonial mesh
 * \return count of components
 * 
 * Simple algorithm that computes the count of components of
 * the Mesh. Isolated vertex is considered as a component
 * Algorithm works as modified BFS
 * 
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
int compute_components(TMesh& m)
{
	typedef TMesh_Traits Mtraits;
	typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
	typedef typename TMesh_Traits::vertex_iterator vertex_iterator;
	typedef typename TMesh_Traits::vv_iterator vv_iterator;

	boost::function_requires<concepts::algorithms::ComputeComponentsConcept<TMesh, TMesh_Traits> >();

	int component_count = 0;

	std::pair<vertex_iterator,vertex_iterator> vertex_pair = Mtraits::get_all_vertices(m);
	std::vector<vertex_descriptor> vertices;

	for (vertex_iterator i = vertex_pair.first;  i != vertex_pair.second; ++i) {
		vertices.push_back(*i);
	}
	std::vector<vertex_descriptor> visited_vertices;

	while (!vertices.empty())
	{
		std::deque<vertex_descriptor> vertex_queue;
		vertex_queue.push_back(*(vertices.begin()));
		visited_vertices.push_back(*(vertices.begin()));

		while (!vertex_queue.empty())
		{
			vertex_descriptor vd = *(vertex_queue.begin());
			vertex_queue.pop_front();

			std::pair<vv_iterator, vv_iterator> adjacent_vertices = Mtraits::get_adjacent_vertices(m, vd);

			for (vv_iterator i = adjacent_vertices.first ; i!= adjacent_vertices.second; ++i)
			{
				vertex_descriptor adjacent_vertex = *i;
				if (std::find(visited_vertices.begin(), visited_vertices.end(), adjacent_vertex) != visited_vertices.end())
					continue;

				vertex_queue.push_front(adjacent_vertex);
				visited_vertices.push_back(adjacent_vertex);
			}
		}

		std::sort(vertices.begin(), vertices.end());
		std::sort(visited_vertices.begin(), visited_vertices.end());

		std::vector<vertex_descriptor> semi_result;

		std::set_difference(vertices.begin(), vertices.end(), visited_vertices.begin(), visited_vertices.end(), std::inserter(semi_result, semi_result.end()));
		vertices = semi_result;
		visited_vertices.clear();

		++component_count;
	}


	return component_count;
}

}//algorithms

}//hmira


#endif /*  __COMPUTE_COMPONENTS_HPP__ */
