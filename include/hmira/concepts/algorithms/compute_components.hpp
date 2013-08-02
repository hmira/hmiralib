#ifdef __COMPUTE_COMPONENTS_HPP__

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

#include <hmira/abstract/mesh_traits.hpp>
#include <hmira/concepts/iterable_vv.h>
#include <hmira/concepts/iterable_vertices.h>

namespace hmira
{

namespace concepts
{


namespace algorithms
{

/*!
 * \defgroup algorithms_concepts The concepts of the algorithms
 * concept check using boost::function_requires
 * 
 * \struct ComputeComponentsConcept
 * \brief concept required for algorithm compute_components
 * 
 * \ingroup algorithms_concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct ComputeComponentsConcept
{
	typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
	typedef typename TMesh_Traits::vertex_iterator vertex_iterator;
	typedef typename TMesh_Traits::vv_iterator vv_iterator;

	void constraints()
	{
		boost::function_requires<hmira::concepts::IterableVVConcept<TMesh, TMesh_Traits> >();
		boost::function_requires<hmira::concepts::IterableVerticesConcept<TMesh, TMesh_Traits> >();
	}
};

}//algorithms

}//concepts

}//hmira
#endif