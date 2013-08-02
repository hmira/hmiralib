#ifdef __MARCHING_CUBES_HPP__

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

#include <hmira/abstract/mesh_traits.hpp>

#include <hmira/concepts/vertex_concept.h>
#include <hmira/concepts/face_concept.h>
#include <hmira/concepts/iterable_faces.h>

namespace hmira
{

namespace concepts
{

/**
 * \brief algorithms
 */
namespace algorithms
{

/*!
 * concept check using boost::function_requires
 * 
 * \struct MarchingCubesConcept
 * \brief concept required for algorithm marching cubes
 * 
 * \ingroup algorithms_concepts
 */
template <class TGrid, class TMesh, class TGrid_Traits, class TMesh_Traits = mesh_traits<TMesh>>
struct MarchingCubesConcept
{
	typedef TGrid						Grid;

	typedef typename TGrid_Traits::Point_descriptor	Point_descriptor;
	typedef typename TGrid_Traits::Point_scalar_type	Point_scalar_type;
	typedef typename TGrid_Traits::Point_properties	Point_properties;
	typedef typename TGrid_Traits::Cube_descriptor		Cube_descriptor;
	typedef typename TGrid_Traits::Cube_iterator		Cube_iterator;

	typedef typename TMesh_Traits::Point			Mesh_Point;
	typedef typename TMesh_Traits::vertex_descriptor	Vertex_descriptor;
	typedef typename TMesh_Traits::face_descriptor		Face_descriptor;
	
	const Grid&	g;
	TMesh&            m;
	Point_descriptor p_d;
	Cube_descriptor cx;
	int i;
	
	void constraints()
	{
		boost::function_requires<hmira::concepts::VertexConcept<TMesh, TMesh_Traits>>();
		boost::function_requires<hmira::concepts::FaceConcept<TMesh, TMesh_Traits>>();
		boost::function_requires<hmira::concepts::IterableFacesConcept<TMesh, TMesh_Traits> >();
		
		Mesh_Point  p( TGrid_Traits::get_coords(g, p_d ));
		TMesh_Traits::get_all_faces(m);
		p_d = TGrid_Traits::get_cube_corner(g, cx, i);
		TGrid_Traits::is_inside( g, p_d );
	}
};


}//algorithms

}//concepts

}//hmira

#endif // __MARCHING_CUBES_HPP__