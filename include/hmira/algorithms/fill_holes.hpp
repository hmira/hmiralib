#ifndef __FILL_HOLES_HPP__
#define __FILL_HOLES_HPP__

#include <hmira/abstract/mesh_traits.hpp>
#include <hmira/concepts/algorithms/fill_holes.hpp>

namespace hmira
{

namespace algorithms
{

/**
 * 
 * \ingroup algorithms
 * 
 * template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>> int fill_holes(TMesh& m)
 * \brief simple algorithm that fills the holes of a mesh
 * \param m a polygonial mesh
 * \return always 0
 * 
 * Simple algorithm that fills holes of a mesh
 * 
 */
template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
int fill_holes(TMesh& m)
{
	boost::function_requires<concepts::algorithms::FillHolesConcept<TMesh, TMesh_Traits> >();

	typedef typename TMesh_Traits::edge_descriptor edge_descriptor;
        typedef typename TMesh_Traits::edge_iterator edge_iterator;

        auto all_edges = TMesh_Traits::get_all_edges(m);

	for (auto i = all_edges.first; i != all_edges.second; ++i)
	{
		edge_descriptor ed = *i;
		if (TMesh_Traits::is_boundary(m, ed)) TMesh_Traits::fill_ring(m, ed);
	}

	return 0;
}

}//algorithms

}//hmira

#endif //__FILL_HOLES_HPP__