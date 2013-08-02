#ifndef ITERABLE_FF_H_
#define ITERABLE_FF_H_

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

namespace hmira
{

namespace concepts
{

/*!
 *  \struct IterableFFConcept
 *  \brief concept required for a operation that iterates through
 *  all faces that are adjacent to a given face
 * 
 *  implemented using boost::function_requires
 *  refines MeshConcept
 * 
 * \ingroup concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct IterableFFConcept
{
	typedef typename TMesh_Traits::face_descriptor face_descriptor;
	typedef typename TMesh_Traits::ff_iterator ff_iterator;

	TMesh m;
	std::pair<ff_iterator,ff_iterator> ffp;
	ff_iterator ffi;
	face_descriptor f;

	void constraints() 
	{
		++ffi;
		f = *ffi;
		ffp = TMesh_Traits::get_adjacent_faces(m, f);
		boost::function_requires< boost::EqualityComparableConcept<ff_iterator> >();
		boost::function_requires< boost::AssignableConcept<ff_iterator> >();
	}
};

}//concepts

}//hmira

#endif /* ITERABLE_FF_H_ */