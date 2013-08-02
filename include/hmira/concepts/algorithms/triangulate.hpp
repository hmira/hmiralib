#ifdef __TRIANGULATE_HPP__

#include <hmira/abstract/mesh_traits.hpp>

namespace hmira
{

namespace concepts
{

namespace algorithms
{

/*!
 * concept check using boost::function_requires
 * 
 * \struct TriangulateConcept
 * \brief concept required for algorithm triangulate
 * 
 * \ingroup algorithms_concepts
 */
template <class TMesh, class TMesh_Traits = mesh_traits<TMesh>>
struct TriangulateConcept
{
	typedef typename TMesh_Traits::face_descriptor face_descriptor;
        typedef typename TMesh_Traits::face_iterator face_iterator;
	typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
	typedef typename TMesh_Traits::fv_iterator fv_iterator;

	face_descriptor fd;
	TMesh& m;
	
	void constraints()
	{
		TMesh_Traits::triangulate_face(m, fd);
	}
};

}//algorithms

}//concepts

}//hmira

#endif //__TRIANGULATE_HPP__