#ifndef EDGE_CONCEPT_H_
#define EDGE_CONCEPT_H_

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
 *  \struct EdgeConcept
 *  \brief concept required for a manipulation with edges
 * 
 *  implemented using boost::function_requires
 * 
 * \ingroup concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct EdgeConcept
{
	typedef typename TMesh_Traits::edge_descriptor edge_descriptor;
	
	TMesh m;
	edge_descriptor e;
	
	void constraints() {
		boost::function_requires< boost::EqualityComparableConcept<edge_descriptor> >();
		boost::function_requires< boost::AssignableConcept<edge_descriptor> >();
	}
};

}//concepts

}//hmira


#endif //EDGE_CONCEPT_H_