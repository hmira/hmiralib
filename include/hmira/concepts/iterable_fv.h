#ifndef ITERABLE_FV_H_
#define ITERABLE_FV_H_

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

namespace hmira
{

namespace concepts
{

/*!
 *  \struct IterableFVConcept
 *  \brief concept required for a operation that iterates through
 *  all vertices that forms a given vertex
 * 
 *  implemented using boost::function_requires
 *  refines MeshConcept
 * 
 * \ingroup concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct IterableFVConcept
{
	typedef typename TMesh_Traits::face_descriptor face_descriptor;
	typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
	typedef typename TMesh_Traits::fv_iterator fv_iterator;

	TMesh m;
	std::pair<fv_iterator,fv_iterator> fvp;
	fv_iterator fvi;
	face_descriptor f;
	vertex_descriptor v;

	void constraints() 
	{
		++fvi;
		v = *fvi;
		fvp = TMesh_Traits::get_surrounding_vertices(m, f);
		boost::function_requires< boost::EqualityComparableConcept<fv_iterator> >();
		boost::function_requires< boost::AssignableConcept<fv_iterator> >();
	}
};

}//concepts

}//hmira

#endif /* ITERABLE_FF_H_ */