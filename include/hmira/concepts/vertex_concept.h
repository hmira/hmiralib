#ifndef VERTEX_CONCEPT_H_
#define VERTEX_CONCEPT_H_

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

namespace hmira
{

namespace concepts
{

/*!
 *  \struct VertexConcept
 *  \brief concept required for a manipulation with vertices
 * 
 *  implemented using boost::function_requires
 * 
 * \ingroup concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct VertexConcept
{
	typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
	typedef typename TMesh_Traits::Point Point;
	
	TMesh m;
	vertex_descriptor v;
	Point p;
	
	void constraints() {
		boost::function_requires< boost::EqualityComparableConcept<vertex_descriptor> >();
		boost::function_requires< boost::AssignableConcept<vertex_descriptor> >();
	}
};


}//concepts

}//hmira


#endif //VERTEX_CONCEPT_H_