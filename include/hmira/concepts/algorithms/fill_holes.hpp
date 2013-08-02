#ifdef __FILL_HOLES_HPP__

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

#include <hmira/abstract/mesh_traits.hpp>

namespace hmira
{

namespace concepts
{

namespace algorithms
{

/*!
 * \defgroup algorithms_concepts
 * concept check using boost::function_requires
 * 
 * \struct FillHolesConcept
 * \brief concept required for algorithm fill_holes
 * 
 * \ingroup algorithms_concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct FillHolesConcept
{
	typedef typename TMesh_Traits::edge_descriptor edge_descriptor;
	typedef typename TMesh_Traits::edge_iterator edge_iterator;

	void constraints()
	{
// 		boost::function_requires<IterableEdgesConcept<TMesh, TMesh_Traits> >();
	}
};

}//algorithms

}//concepts

}//hmira

#endif //__FILL_HOLES_HPP__