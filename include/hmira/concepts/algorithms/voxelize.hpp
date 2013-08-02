#ifdef __VOXELIZE_HPP__

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

#include <hmira/abstract/mesh_traits.hpp>

#include <hmira/concepts/iterable_faces.h>

namespace hmira
{

namespace concepts
{

namespace algorithms
{

/*!
 * \defgroup algorithms_concepts The concepts of the Algorithms
 * concept check using boost::function_requires
 * 
 * \struct TriangulateConcept
 * \brief concept required for algorithm voxelize
 * 
 * \ingroup algorithms_concepts
 */
template <
	typename TGrid,
	typename TMesh,
	typename TGrid_Traits,// = grid_traits<TGrid>, 
	typename TMesh_Traits = mesh_traits<TMesh>
	>
struct VoxelizeConcept
{
	typedef TGrid						Grid;

	typedef typename TGrid_Traits::Point_descriptor	Point_descriptor;
	typedef typename TGrid_Traits::Coordinates_type	Coordinates_descriptor;
	typedef typename TGrid_Traits::Point_scalar_type	Point_scalar_type;
	typedef typename TGrid_Traits::Point_properties	Point_properties;
	typedef typename TGrid_Traits::Cube_descriptor		Cube_descriptor;
	typedef typename TGrid_Traits::Cube_iterator		Cube_iterator;
	typedef typename TGrid_Traits::Vector_unit_type	Grid_unit_type;

	typedef typename TMesh_Traits::Point			Mesh_Point;
	typedef typename TMesh_Traits::vertex_descriptor	Vertex_descriptor;
	typedef typename TMesh_Traits::face_descriptor		Face_descriptor;
	
	Grid& g;
	Grid_unit_type elem;
	TMesh& m;
	Face_descriptor f;
	Mesh_Point mp;
	Point_descriptor pd;
	Coordinates_descriptor cd;
	Vertex_descriptor vd;
	
	void constraints()
	{
		boost::function_requires<hmira::concepts::VertexConcept<TMesh, TMesh_Traits>>();
		boost::function_requires<hmira::concepts::FaceConcept<TMesh, TMesh_Traits>>();
		boost::function_requires<hmira::concepts::IterableFacesConcept<TMesh, TMesh_Traits> >();
		
		elem = TGrid_Traits::get_x_element_size(g);
		elem = TGrid_Traits::get_y_element_size(g);
		elem = TGrid_Traits::get_z_element_size(g);
		mp = TMesh_Traits::get_coordinates(m, vd);
		
		Coordinates_descriptor v(0.f, 0.f, 0.f);
		Cube_descriptor pt = TGrid_Traits::get_cube_from_coords(g, v);
		TGrid_Traits::fill_cube_by_value(g, pt, -1);
		TGrid_Traits::is_cube_inside(g, pt);
	}
};

}//algorithms

}//concepts

}//hmira

#endif //__VOXELIZE_HPP__