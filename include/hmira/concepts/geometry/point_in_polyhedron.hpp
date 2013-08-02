#ifdef __POINT_IN_POLYHEDRON_HPP__

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

#include <hmira/abstract/mesh_traits.hpp>
#include <hmira/concepts/iterable_faces.h>
#include <hmira/concepts/face_concept.h>

namespace hmira
{

/**
 * 
 * 
 * \ingroup concepts
 */
namespace concepts
{

/**
 * 
 * \brief A namespace algorithms
 * \ingroup algorithms
 */

namespace geometry
{

/*!
 * \defgroup algorithms_concepts
 * concept check using boost::function_requires
 * 
 * \struct PointInPolyhedronConcept
 * \brief concept required for algorithm that determines whether the point
 * is either in or out or even on the surface
 * 
 * \ingroup algorithms_concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct PointInPolyhedronConcept
{
	typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;

	void constraints()
	{
		boost::function_requires<hmira::concepts::IterableFacesConcept<TMesh_Traits>>();
	}
};

}//geometry

}//concepts

}//hmira
#endif
