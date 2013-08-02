#ifndef __FLIP_NORMALS_HPP__
#define __FLIP_NORMALS_HPP__

#include <hmira/abstract/mesh_traits.hpp>
#include <hmira/concepts/algorithms/flip_normals.hpp>

namespace hmira
{

namespace algorithms
{

/**
 * 
 * \ingroup algorithms
 * 
 * \fn template <class TMesh, class TMesh_Traits> int flip_normals(TMesh& m)
 * \brief flips normals
 * \param m a polygonial mesh
 * \return always 0
 * 
 * Algorithm that iterates through all faces
 * and flips its normal
 * 
 */
template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
int flip_normals(TMesh& m)
{
	typedef typename TMesh_Traits::face_descriptor face_descriptor;
	typedef typename TMesh_Traits::face_iterator face_iterator;

	boost::function_requires<concepts::algorithms::FlipNormalsConcept<TMesh, TMesh_Traits> >();
	
	auto all_faces = TMesh_Traits::get_all_faces(m);

	for (auto i = all_faces.first; i != all_faces.second; ++i)
	{
		face_descriptor fd = *i;
		TMesh_Traits::flip_face_normal(m, fd);
	}

	return 0;
}

} //algorithms

} //hmira

#endif //__FLIP_NORMALS_HPP__