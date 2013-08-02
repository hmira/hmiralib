#ifndef FACE_CONCEPT_H_
#define FACE_CONCEPT_H_

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

/**
 * \brief A top namespace of the library
 * 
 */
namespace hmira
{
/**
 * \brief A namespace with concepts that check
 * whether the traits are build correctly for an algorithm
 * 
 */

namespace concepts
{
/*!
 *  \struct FaceConcept
 *  \brief concept required for a manipulation with faces
 * 
 *  implemented using boost::function_requires
 * 
 * \ingroup concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct FaceConcept
{
	typedef typename TMesh_Traits::face_descriptor face_descriptor;
	
	TMesh m;
	face_descriptor f;
	
	void constraints() {
		boost::function_requires< boost::EqualityComparableConcept<face_descriptor> >();
	}
};

}//concepts

}//hmira


#endif //FACE_CONCEPT_H_
