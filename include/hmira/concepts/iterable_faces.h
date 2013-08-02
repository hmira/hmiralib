/*
 * iterable_vertices.h
 *
 *  Created on: Jul 24, 2012
 *      Author: hmirap
 */

#ifndef ITERABLE_FACES_H_
#define ITERABLE_FACES_H_

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

namespace hmira
{

namespace concepts
{

/*!
 *  \struct IterableFacesConcept
 *  \brief concept required for a operation that iterates through all faces
 * 
 *  implemented using boost::function_requires
 * 
 * \ingroup concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct IterableFacesConcept
{
	typedef typename TMesh_Traits::face_descriptor face_descriptor;
	typedef typename TMesh_Traits::face_iterator face_iterator;

	TMesh m;
	std::pair<face_iterator,face_iterator> fp;
	face_descriptor f;
	face_iterator fi;
	
	void constraints() {
		
		f = *fi;
		++fi;
		
		fp = TMesh_Traits::get_all_faces(m);
	}
};

}//concepts

}//hmira

#endif /* ITERABLE_FACES_H_ */
