#ifndef ITERABLE_VV_H_
#define ITERABLE_VV_H_

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

namespace hmira
{

namespace concepts
{

/*!
 *  \struct IterableVVConcept
 *  \brief concept required for a operation that iterates through
 *  all vertices that are adjacent to a given vertex
 * 
 *  implemented using boost::function_requires
 *  refines MeshConcept
 * 
 * \ingroup concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct IterableVVConcept
{
	typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
	typedef typename TMesh_Traits::vv_iterator vv_iterator;

	TMesh m;
	std::pair<vv_iterator,vv_iterator> vvp;
	vv_iterator vvi;
	vertex_descriptor v;

	void constraints() 
	{
		++vvi;
		v = *vvi;
		vvp = TMesh_Traits::get_adjacent_vertices(m, v);
		boost::function_requires< boost::EqualityComparableConcept<vv_iterator> >();
		boost::function_requires< boost::AssignableConcept<vv_iterator> >();
	}
};

}//concepts

}//hmira

#endif /* ITERABLE_VV_H_ */