#ifdef __FLIP_NORMALS_HPP__

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

#include <hmira/abstract/mesh_traits.hpp>

namespace hmira
{

/**
 * 
 * \brief concepts
 */
namespace concepts
{

/**
 * 
 * \brief algorithms
 */

namespace algorithms
{

/*!
 * concept check using boost::function_requires
 * 
 * \struct FlipNormalsConcept
 * \brief concept required for algorithm flip_face_normal
 * 
 * \ingroup algorithms_concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct FlipNormalsConcept
{
	typedef typename TMesh_Traits::face_descriptor face_descriptor;
	typedef typename TMesh_Traits::face_iterator face_iterator;

	TMesh& m;
	face_descriptor f;
	
	void constraints()
	{
		TMesh_Traits::get_all_faces(m);
		TMesh_Traits::flip_face_normal(m, f);
// 		boost::function_requires<IterableFacesConcept<TMesh, TMesh_Traits> >();
	}
};

}//algorithms

}//concepts

}//hmira

#endif //__FLIP_NORMALS_HPP__