#ifndef __TRIANGULATE_HPP__
#define __TRIANGULATE_HPP__

#include <hmira/abstract/mesh_traits.hpp>
#include <hmira/concepts/algorithms/triangulate.hpp>

namespace hmira
{

namespace algorithms
{

/**
 * 
 * \ingroup algorithms
 * 
 * template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>> int triangulate(TMesh& m)
 * \brief flips normals
 * \param m a polygonial mesh
 * \return always 0
 * 
 * Simple algorithm that iterates through all faces
 * and triangulate it
 * 
 */
template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
int triangulate(TMesh& m)
{
	typedef typename TMesh_Traits::face_descriptor face_descriptor;
        typedef typename TMesh_Traits::face_iterator face_iterator;
	typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
	typedef typename TMesh_Traits::fv_iterator fv_iterator;

	boost::function_requires<hmira::concepts::algorithms::TriangulateConcept<TMesh, TMesh_Traits> >();
	
        auto all_faces = TMesh_Traits::get_all_faces(m);

	for (auto i = all_faces.first; i != all_faces.second; ++i)
	{
		face_descriptor fd = *i;
		TMesh_Traits::triangulate_face(m, fd);
	}

	return 0;
}

}//algorithms

}//hmira

#endif //__TRIANGULATE_HPP__