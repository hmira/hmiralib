#ifndef __VOXELIZE_HPP__
#define __VOXELIZE_HPP__

#include <limits>
#include <hmira/geometry/ray_face_intersection.hpp>

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

#include <hmira/abstract/mesh_traits.hpp>
#include <hmira/concepts/algorithms/voxelize.hpp>

namespace hmira
{
/**
 * \brief A namespace algorithms
 * 
 */
namespace algorithms
{

/**
 * \defgroup	algorithms Algorithms
 * 
 * \ingroup algorithms
 * 
 * \class Voxelize
 * \brief converts the mesh into the voxel representation
 * 
 * voxelization algorithm
 * 
 * the algorithm first process a edges, then fills the faces
 * and finally fills the object line-by-line
 * 
 * <b>the triangle mesh is required. If the non-triangle mesh
 * is used, the algorithm does not work correctly</b>
 * 
 */
template <
	typename TGrid,
	typename TMesh,
	typename TGrid_Traits,// = grid_traits<TGrid>, 
	typename TMesh_Traits = mesh_traits<TMesh>
	>
class Voxelize
{

public:
	
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

	Grid&	grid_;
	Grid_unit_type elem_x_, elem_y_, elem_z_;
	TMesh&            mesh_;
	
	~Voxelize(){}
	Voxelize( Grid& _grid, TMesh& _mesh)
		: grid_( _grid ),
		mesh_( _mesh )
	{
		boost::function_requires< concepts::algorithms::VoxelizeConcept<TGrid, TMesh, TGrid_Traits, TMesh_Traits> >();

		elem_x_ = TGrid_Traits::get_x_element_size(grid_);
		elem_y_ = TGrid_Traits::get_y_element_size(grid_);
		elem_z_ = TGrid_Traits::get_z_element_size(grid_);
	}
	
	/**
	 * the method that process all faces and build the voxelized skeleton
	 * 
	 */
	void process_edges()
	{
		auto all_faces = TMesh_Traits::get_all_faces(mesh_);
		for (auto f = all_faces.first; f != all_faces.second; ++f)
		{
			auto minX = std::numeric_limits< float >::max();
			auto maxX = std::numeric_limits< float >::min();
			auto minY = std::numeric_limits< float >::max();
			auto maxY = std::numeric_limits< float >::min();
			auto minZ = std::numeric_limits< float >::max();
			auto maxZ = std::numeric_limits< float >::min();
			
			auto f_vertex_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
			for (auto fv_i = f_vertex_pair.first; fv_i != f_vertex_pair.second; ++fv_i)
			{
				auto vertex = *fv_i;
				auto coords = TMesh_Traits::get_coordinates(mesh_, vertex);
				if (minX > coords[0]) minX = coords[0];
				if (maxX < coords[0]) maxX = coords[0];
				if (minY > coords[1]) minY = coords[1];
				if (maxY < coords[1]) maxY = coords[1];
				if (minZ > coords[2]) minZ = coords[2];
				if (maxZ < coords[2]) maxZ = coords[2];
			}
			
			auto dx = maxX - minX;
			auto dy = maxY - minY;
			auto dz = maxZ - minZ;
			
			rasterize_triangle(*f);
		}
	}

	/**
	 * the method that process all faces and build the voxelized surface
	 * according the skeleton build by the method <i>process_edges()</i>
	 * 
	 * fills the each face according the pre-computed slope the plane that
	 * contains the face
	 */
	void process_faces()
	{
		auto all_faces = TMesh_Traits::get_all_faces(mesh_);
		for (auto f = all_faces.first; f != all_faces.second; ++f)
		{
			auto minX = std::numeric_limits< float >::max();
			auto maxX = std::numeric_limits< float >::min();
			auto minY = std::numeric_limits< float >::max();
			auto maxY = std::numeric_limits< float >::min();
			auto minZ = std::numeric_limits< float >::max();
			auto maxZ = std::numeric_limits< float >::min();
			
			auto f_vertex_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
			for (auto fv_i = f_vertex_pair.first; fv_i != f_vertex_pair.second; ++fv_i)
			{
				auto vertex = *fv_i;
				auto coords = TMesh_Traits::get_coordinates(mesh_, vertex);
				if (minX > coords[0]) minX = coords[0];
				if (maxX < coords[0]) maxX = coords[0];
				if (minY > coords[1]) minY = coords[1];
				if (maxY < coords[1]) maxY = coords[1];
				if (minZ > coords[2]) minZ = coords[2];
				if (maxZ < coords[2]) maxZ = coords[2];
			}
			
			auto dx = maxX - minX;
			auto dy = maxY - minY;
			auto dz = maxZ - minZ;
			
			rasterize_triangle(*f);
			
			if (fabs(dz) < fabs(dx) && fabs(dz) < fabs(dy))
				fill_triangle_xy(*f);
			else if (fabs(dy) < fabs(dx) && fabs(dy) < fabs(dz))
				fill_triangle_xz(*f);
			else
				fill_triangle_yz(*f);
		}
	}

	/**
	 * the method that process all faces and build the voxelized surface
	 * according the skeleton build by the method <i>process_edges()</i>
	 * 
	 * finally, it runs the floodfill over the voxelized structure
	 */
	void process_volume()
	{
		this->process_faces();
		this->fill_space();
	}
	
	void process_cube(Cube_descriptor cube)
	{
		Point_descriptor	corner[8];
		Vertex_descriptor	samples[12];
		unsigned int		cubetype( 0 );
		unsigned int		i;

		for ( i=0; i<8; ++i )
			corner[i] = TGrid_Traits::get_cube_corner(grid_, cube, i );

	}

	/**
	 * voxelize all edges of the triangle
	 * 
	 * \param	f the triangle to be voxelized
	 */
	bool rasterize_triangle(Face_descriptor f)
	{
		auto f_edge_pair = TMesh_Traits::get_surrounding_edges(mesh_, f);
		for (auto fe_i = f_edge_pair.first; fe_i != f_edge_pair.second; ++fe_i)
		{
			auto v_pair = TMesh_Traits::get_edge_vertices(mesh_, *fe_i);
			auto v0 = v_pair.first;
			auto v1 = v_pair.second;
			
			rasterize_line(v0, v1);
			
		}
		
		
	}

	/**
	 * \brief rasterizes the edge given by the pair of the vertices
	 * 
	 * first, pre-compute a slope of the line and then voxelize the line
	 * 
	 * \param	v0 first vertex
	 * \param	v1 second vertex
	 */
	bool rasterize_line(Vertex_descriptor v0, Vertex_descriptor v1)
	{
		auto p0 = TMesh_Traits::get_coordinates(mesh_, v0);
		auto p1 = TMesh_Traits::get_coordinates(mesh_, v1);
		
		auto p0x = p0[0];
		auto p0y = p0[1];
		auto p0z = p0[2];
		
		auto p1x = p1[0];
		auto p1y = p1[1];
		auto p1z = p1[2];
		
		auto dx = p0x - p1x;
		auto dy = p0y - p1y;
		auto dz = p0z - p1z;
		
		if (fabs(dx) > fabs(dy) && fabs(dx) > fabs(dz))
		{
			rasterize_line_x(p0, p1);
		}
		
		else if (fabs(dy) > fabs(dx) && fabs(dy) > fabs(dz))
		{
			rasterize_line_y(p0, p1);
		}
		
		else if (fabs(dz) > fabs(dx) && fabs(dz) > fabs(dy))
		{
			rasterize_line_z(p0, p1);
		}
	}
	
	bool rasterize_line_x(Mesh_Point &p0, Mesh_Point &p1)
	{	
		auto p0x = p0[0];
		auto p0y = p0[1];
		auto p0z = p0[2];
		
		auto p1x = p1[0];
		auto p1y = p1[1];
		auto p1z = p1[2];
		
		int ninja_constant = (p0x < p1x) ? 1 : -1;
		
		auto dx = p0x - p1x;
		auto dy = p0y - p1y;
		auto dz = p0z - p1z;
		
		auto dyx = dy / dx;
		auto dzx = dz / dx;
		
		auto x = p0x;
		auto y = p0y;
		auto z = p0z;
		
		while ( ( x < p1x - elem_x_ ) || (x > p1x + elem_x_ ))
		{
			Coordinates_descriptor v(x,y,z);
			auto pt = TGrid_Traits::get_cube_from_coords(grid_, v);
			TGrid_Traits::fill_cube_by_value(grid_, pt, -1);
			
			x += ninja_constant * elem_x_;
			y += ninja_constant * dyx * elem_x_;
			z += ninja_constant * dzx * elem_x_;
		}
		
		Coordinates_descriptor v(x,y,z);
		auto pt = TGrid_Traits::get_cube_from_coords(grid_, v);
		TGrid_Traits::fill_cube_by_value(grid_, pt, -1);
	}
	bool rasterize_line_y(Mesh_Point &p0, Mesh_Point &p1)
	{		
		auto p0x = p0[0];
		auto p0y = p0[1];
		auto p0z = p0[2];
		
		auto p1x = p1[0];
		auto p1y = p1[1];
		auto p1z = p1[2];
		
		int ninja_constant = (p0y < p1y) ? 1 : -1;
		
		auto dx = p0x - p1x;
		auto dy = p0y - p1y;
		auto dz = p0z - p1z;
		
		auto dxy = dx / dy;
		auto dzy = dz / dy;
		
		auto x = p0x;
		auto y = p0y;
		auto z = p0z;
		
		while ( (y < p1y - elem_y_) || (y > p1y + elem_y_ ))
		{
			Coordinates_descriptor v(x,y,z);
			auto pt = TGrid_Traits::get_cube_from_coords(grid_, v);
			TGrid_Traits::fill_cube_by_value(grid_, pt, -1);
			
			x += ninja_constant * dxy * elem_y_;
			y += ninja_constant * elem_y_;
			z += ninja_constant * dzy * elem_y_;
		}
		
		Coordinates_descriptor v(x,y,z);
		auto pt = TGrid_Traits::get_cube_from_coords(grid_, v);
		TGrid_Traits::fill_cube_by_value(grid_, pt, -1);
	}
	bool rasterize_line_z(Mesh_Point &p0, Mesh_Point &p1)
	{
		auto p0x = p0[0];
		auto p0y = p0[1];
		auto p0z = p0[2];
		
		auto p1x = p1[0];
		auto p1y = p1[1];
		auto p1z = p1[2];
		
		int ninja_constant = (p0z < p1z) ? 1 : -1;
		
		auto dx = p0x - p1x;
		auto dy = p0y - p1y;
		auto dz = p0z - p1z;
		
		auto dxz = dx / dz;
		auto dyz = dy / dz;
		
		auto x = p0x;
		auto y = p0y;
		auto z = p0z;
		
		while ( (z < p1z - elem_z_) || (z > p1z + elem_z_))
		{
			Coordinates_descriptor v(x,y,z);
			auto pt = TGrid_Traits::get_cube_from_coords(grid_, v);
			TGrid_Traits::fill_cube_by_value(grid_, pt, -1);
			
			x += ninja_constant * dxz * elem_z_;
			y += ninja_constant * dyz * elem_z_;
			z += ninja_constant * elem_z_;
		}
		
		Coordinates_descriptor v(x,y,z);
		auto pt = TGrid_Traits::get_cube_from_coords(grid_, v);
		TGrid_Traits::fill_cube_by_value(grid_, pt, -1);
	}

	/**
	 * fills the cube according the <b>inside status</b>
	 * 
	 * \param	x x-coordinate of the cube
	 * \param	y y-coordinate of the cube
	 * \param	z z-coordinate of the cube
	 * \param	inside inside status
	 * \param	cubes the set of cubes
	 */
	bool fill_cube_one(
		float x,
		float y,
		float z,
		int& inside,
		std::vector<Cube_descriptor>& cubes
	)
	{
		if ((x > 0) && (y > 0) && (z > 0))
		{
			Coordinates_descriptor v(x,y,z);
			auto pt = TGrid_Traits::get_cube_from_coords(grid_, v);
			if (TGrid_Traits::is_cube_inside(grid_,pt))
			{
				if (inside == 0) //in empty
				{
					inside = 2; //awaiting fill
				}
				else if (inside == 1) //awaiting empty
				{
					inside = 1; //empty space
				}
				else if (/*inside == 2 || */inside == 3)
				{
					inside = 1;
				}
			}
			else
			{
				if (inside == 1)
				{
					inside = 0;
				}
				else if (inside == 2) //awaiting fill
				{
					cubes.push_back(pt);
					inside = 3; //inside object
				}
				else if (inside == 3) //inside object
				{
					cubes.push_back(pt);
				}
			}
		}
	}

	/**
	 * fills the cube according the <b>inside status</b>
	 * cubes are placed into the set of cubes <b>cubes</b>
	 * that is later processed
	 * 
	 * \param	cube the given cube
	 * \param	inside inside status
	 * \param	cubes the set of cubes to be filled
	 */
	bool fill_cube_one(
		Cube_descriptor cube,
		int& inside,
		std::vector<Cube_descriptor>& cubes)
	{
		if (TGrid_Traits::is_cube_inside(grid_, cube))
		{
			if (inside == 0) //in empty
			{
				inside = 2; //awaiting fill
			}
			else if (inside == 1) //awaiting empty
			{
				inside = 1; //empty space
			}
			else if (/*inside == 2 ||*/ inside == 3)
			{
				inside = 1;
			}
		}
		else
		{
			if (inside == 1)
			{
				inside = 0;
			}
			else if (inside == 2) //awaiting fill
			{
				cubes.push_back(cube);
				inside = 3; //inside object
			}
			else if (inside == 3) //inside object
			{
				cubes.push_back(cube);
			}
		}
	}

	/**
	 * \param	inside the inside status
	 * \param	cubes	the set of cubes to be filled
	 * 
	 * \brief the purpose of this function is to avoid to fill the rows that
	 * crosses the odd count of faces
	 * 
	 * if the inside status does not indicate that the floodfill process is outside
	 * the row can not be filled
	 * 
	 */
	bool apply_fill(int inside, std::vector<Cube_descriptor>& cubes)
	{
		if (inside == 3)
		{
			cubes.clear();
			return false;
		}
		for (auto cube : cubes)
		{
			TGrid_Traits::fill_cube_by_value(grid_, cube, -1);			
		}
		cubes.clear();
		return true;
	}
	
	bool fill_triangle_xy(Face_descriptor f)
	{
		auto minX = std::numeric_limits< float >::max();
		auto maxX = std::numeric_limits< float >::min();
		auto minY = std::numeric_limits< float >::max();
		auto maxY = std::numeric_limits< float >::min();
		auto minZ = std::numeric_limits< float >::max();
		auto maxZ = std::numeric_limits< float >::min();
		
		Coordinates_descriptor coords;
		
		auto f_vertex_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
		for (auto fv_i = f_vertex_pair.first; fv_i != f_vertex_pair.second; ++fv_i)
		{
			auto vertex = *fv_i;
			coords = TMesh_Traits::get_coordinates(mesh_, vertex);
			if (minX > coords[0]) minX = coords[0];
			if (maxX < coords[0]) maxX = coords[0];
			if (minY > coords[1]) minY = coords[1];
			if (maxY < coords[1]) maxY = coords[1];
			if (minZ > coords[2]) minZ = coords[2];
			if (maxZ < coords[2]) maxZ = coords[2];
		}
		
		auto norm = TMesh_Traits::get_face_normal(mesh_, f); //cross(v0, v1);
		
		/*parameters*/
		auto a = norm[0];
		auto b = norm[1];
		auto c = norm[2];
		auto d = -1 * (a * coords[0] + b * coords[1] + c * coords[2]);
		
		auto dx = maxX - minX;
		auto dy = maxY - minY;

		auto dyx = dy / dx;
		auto dxy = dx / dy;
		
		auto x = minX;
		auto y = minY;
		
		bool parallel_to_axis = (c < 0.00001 && c > -0.00001);
		
		auto rcp_c = 1.f / c;
		auto z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
		int inside = 0;
		std::vector<Cube_descriptor> cubes;
		if (dx > dy)
		{
			while (( y < maxY - dyx * elem_x_) || (y > maxY + dyx * elem_x_))
			{
				inside = 0;
				x = minX;
				z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
				
				while (( x < maxX - elem_x_ ) || (x > maxX + elem_x_ ))
				{
					fill_cube_one(x,y,z,inside,cubes);
					x += elem_x_;
					z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
				}
				fill_cube_one(x,y,z,inside, cubes);
				apply_fill(inside, cubes);
				y += dyx * elem_x_;
			}
			x = minX;
			z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
			
			inside = 0;
			while (( x < maxX - elem_x_ ) || (x > maxX + elem_x_ ))
			{
				fill_cube_one(x,y,z,inside,cubes);
				x += elem_x_;
				z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
			}
			fill_cube_one(x,y,z,inside,cubes);
			apply_fill(inside, cubes);

		}
		else
		{
			while (( x < maxX - dxy * elem_y_) || (x > maxX + dxy * elem_y_))
			{
				y = minY;
				z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
				
				inside = 0;
				while (( y < maxY - elem_y_ ) || (y > maxY + elem_y_ ))
				{
					fill_cube_one(x,y,z,inside,cubes);
					y += elem_y_;
					z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
				}
				fill_cube_one(x,y,z,inside,cubes);
				apply_fill(inside, cubes);
				x += dxy * elem_y_;
			}
			
			y = minY;
			z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
			
			inside = 0;
			while (( y < maxY - elem_y_ ) || (y > maxY + elem_y_ ))
			{
				fill_cube_one(x,y,z,inside,cubes);
				y += elem_y_;
				z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
			}
			fill_cube_one(x,y,z,inside,cubes);
			apply_fill(inside, cubes);

		}
	}
	
	bool fill_triangle_xz(Face_descriptor f)
	{
		auto minX = std::numeric_limits< float >::max();
		auto maxX = std::numeric_limits< float >::min();
		auto minY = std::numeric_limits< float >::max();
		auto maxY = std::numeric_limits< float >::min();
		auto minZ = std::numeric_limits< float >::max();
		auto maxZ = std::numeric_limits< float >::min();
		
		Coordinates_descriptor coords;
		
		auto f_vertex_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
		for (auto fv_i = f_vertex_pair.first; fv_i != f_vertex_pair.second; ++fv_i)
		{
			auto vertex = *fv_i;
			coords = TMesh_Traits::get_coordinates(mesh_, vertex);
			if (minX > coords[0]) minX = coords[0];
			if (maxX < coords[0]) maxX = coords[0];
			if (minY > coords[1]) minY = coords[1];
			if (maxY < coords[1]) maxY = coords[1];
			if (minZ > coords[2]) minZ = coords[2];
			if (maxZ < coords[2]) maxZ = coords[2];
		}
		
		auto norm = TMesh_Traits::get_face_normal(mesh_, f); //cross(v0, v1);
		
		/*parameters*/
		auto a = norm[0];
		auto b = norm[1];
		auto c = norm[2];
		auto d = -1 * (a * coords[0] + b * coords[1] + c * coords[2]);
		
		auto dx = maxX - minX;
		auto dz = maxZ - minZ;

		auto dzx = dz / dx;
		auto dxz = dx / dz;
		
		auto x = minX;
		auto z = minZ;
		
		bool parallel_to_axis = (b < 0.00001 && b > -0.00001);
		auto rcp_b = 1.f / b;
		auto y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
		
		int inside = 0;
		std::vector<Cube_descriptor> cubes;
		
		if (dx > dz)
		{
			while (( z < maxZ - dzx * elem_x_) || (z > maxZ + dzx * elem_x_))
			{
				inside = 0;
				x = minX;
				y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
				
				while (( x < maxX - elem_x_ ) || (x > maxX + elem_x_ ))
				{
					fill_cube_one(x,y,z,inside,cubes);
					x += elem_x_;
					y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
				}
				fill_cube_one(x,y,z,inside,cubes);
				apply_fill(inside, cubes);
				z += dzx * elem_x_;
			}
			/*LAST ROW*/
			x = minX;
			y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
			inside = 0;
			
			while (( x < maxX - elem_x_ ) || (x > maxX + elem_x_ ))
			{
				fill_cube_one(x,y,z,inside,cubes);
				x += elem_x_;
				y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
			}
			fill_cube_one(x,y,z,inside,cubes);
			apply_fill(inside, cubes);

		}
		else
		{
			while (( x < maxX - dxz * elem_z_) || (x > maxX + dxz * elem_z_))
			{
				inside = 0;
				z = minZ;
				y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
				
				while (( z < maxZ - elem_z_ ) || (z > maxZ + elem_z_ ))
				{
					fill_cube_one(x,y,z,inside,cubes);
					z += elem_z_;
					y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
				}
				fill_cube_one(x,y,z,inside,cubes);
				apply_fill(inside, cubes);
				x += dxz * elem_z_;
			}
			
			z = minZ;
			y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
			inside = 0;
			
			while (( z < maxZ - elem_z_ ) || (z > maxZ + elem_z_ ))
			{
				fill_cube_one(x,y,z,inside,cubes);
				z += elem_z_;
				y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
			}
			fill_cube_one(x,y,z,inside,cubes);
			apply_fill(inside, cubes);
		}
	}
	
	bool fill_triangle_yz(Face_descriptor f)
	{
		auto minX = std::numeric_limits< float >::max();
		auto maxX = std::numeric_limits< float >::min();
		auto minY = std::numeric_limits< float >::max();
		auto maxY = std::numeric_limits< float >::min();
		auto minZ = std::numeric_limits< float >::max();
		auto maxZ = std::numeric_limits< float >::min();
		
		Coordinates_descriptor coords;
		
		auto f_vertex_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
		for (auto fv_i = f_vertex_pair.first; fv_i != f_vertex_pair.second; ++fv_i)
		{
			auto vertex = *fv_i;
			coords = TMesh_Traits::get_coordinates(mesh_, vertex);
			if (minX > coords[0]) minX = coords[0];
			if (maxX < coords[0]) maxX = coords[0];
			if (minY > coords[1]) minY = coords[1];
			if (maxY < coords[1]) maxY = coords[1];
			if (minZ > coords[2]) minZ = coords[2];
			if (maxZ < coords[2]) maxZ = coords[2];
		}
		
		auto norm = TMesh_Traits::get_face_normal(mesh_, f); //cross(v0, v1);
		
		/*parameters*/
		auto a = norm[0];
		auto b = norm[1];
		auto c = norm[2];
		auto d = -1 * (a * coords[0] + b * coords[1] + c * coords[2]);
		
		auto dz = maxZ - minZ;
		auto dy = maxY - minY;

		auto dyz = dy / dz;
		auto dzy = dz / dy;
		
		auto z = minZ;
		auto y = minY;
		bool parallel_to_axis = (a < 0.00001 && a > -0.00001);
		auto rcp_a = 1.f / a;
		auto x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
		
		int inside = 0;
		std::vector<Cube_descriptor> cubes;
		
		if (dz > dy)
		{
			while (( y < maxY - dyz * elem_z_) || (y > maxY + dyz * elem_z_))
			{
				inside = 0;
				z = minZ;
				x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
				
				while (( z < maxZ - elem_z_ ) || ( z > maxZ + elem_z_ ))
				{
					fill_cube_one(x,y,z,inside,cubes);
					z += elem_z_;
					x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
				}
				fill_cube_one(x,y,z,inside,cubes);
				apply_fill(inside, cubes);
				y += dyz * elem_z_;
			}
			z = minZ;
			x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
			inside = 0;
			
			while (( z < maxZ - elem_z_ ) || ( z > maxZ + elem_z_ ))
			{
				fill_cube_one(x,y,z,inside,cubes);
				z += elem_z_;
				x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
			}
			fill_cube_one(x,y,z,inside,cubes);
			apply_fill(inside, cubes);
		}
		else
		{
			while (( z < maxZ - dzy * elem_y_) || (z > maxZ + dzy * elem_y_))
			{
				inside = 0;
				y = minY;
				x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
				
				while (( y < maxY - elem_y_ ) || ( y > maxY + elem_y_ ))
				{
					fill_cube_one(x,y,z,inside,cubes);
					y += elem_y_;
					x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
				}
				fill_cube_one(x,y,z,inside,cubes);
				apply_fill(inside, cubes);
				z += dzy * elem_y_;
			}
			y = minY;
			x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
			inside = 0;
			
			while (( y < maxY - elem_y_ ) || ( y > maxY + elem_y_ ))
			{
				fill_cube_one(x,y,z,inside,cubes);
				y += elem_y_;
				x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
			}
			fill_cube_one(x,y,z,inside,cubes);
			apply_fill(inside, cubes);

		}
	}

	/**
	 * floodfill method
	 */
	bool
	fill_space()
	{
		int x = 0;
		int y = 0;
		int z = 0;
		
		auto maxX = TGrid_Traits::get_x_resolution(grid_);
		auto maxY = TGrid_Traits::get_y_resolution(grid_);
		auto maxZ = TGrid_Traits::get_z_resolution(grid_);
		
		/**
		 * inside status
		 */
		int inside = 0;
		int i = 0;
		std::vector<Cube_descriptor> cubes;
		for (auto cube : grid_)
		{
			if ( i % maxX == 0)
			{
				apply_fill(inside, cubes);
				inside = 0;
			}
			fill_cube_one(cube, inside, cubes);
			i++;
		}
	}
	
	
	bool crosses_face(Point_descriptor corner0, Point_descriptor corner1, Face_descriptor f)
	{
		auto c0 = TGrid_Traits::get_coords(grid_, corner0);
		auto c1 = TGrid_Traits::get_coords(grid_, corner1);
		auto fv_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
		Vertex_descriptor v[3];
		
		auto fv = fv_pair.first;
		for(int i = 0; i < 3; i++ )
		{
			v[i] = *fv;
			++fv;
		}
		
		auto a = (Coordinates_descriptor)TMesh_Traits::get_coordinates(mesh_, v[0]);
		auto b = (Coordinates_descriptor)TMesh_Traits::get_coordinates(mesh_, v[1]);
		auto c = (Coordinates_descriptor)TMesh_Traits::get_coordinates(mesh_, v[2]);
		
		auto c_direction = c1 - c0;
		
		Coordinates_descriptor intersection;
		float distance = 0.f;
		
		auto result = line_face_intersection(
		c0,
		c_direction,
		a,
		b,
		c,
		intersection,
		distance);
		
		if ((c0 - c1).length() < distance)
			result = false;
		
		
		if (result)
		{
			auto normal = TMesh_Traits::get_face_normal(mesh_, f);
			if (dot(normal, c_direction) > 0)
			{
				TGrid_Traits::set_scalar_value( grid_, corner0, 0.5);
				if (TGrid_Traits::scalar_value( grid_, corner1) == -0.5)
				{
					TGrid_Traits::set_scalar_value( grid_, corner1, -0.5);
				}
				else
				{
					TGrid_Traits::set_scalar_value( grid_, corner1, -0.5);
				}
			}
			else
			{
				TGrid_Traits::set_scalar_value( grid_, corner1, 0.5);
				if (TGrid_Traits::scalar_value( grid_, corner0) == -0.5)
				{
					TGrid_Traits::set_scalar_value( grid_, corner0, -0.5);
				}
				else
				{
					TGrid_Traits::set_scalar_value( grid_, corner0, -0.5);
				}
			}
		}
		return result;
	}

	/**
	 * \brief set cube as filled
	 * 
	 * \param	_cidx cube to be filled
	 */
	void fill_cube(unsigned int _cidx)
	{
		auto corner0 = TGrid_Traits::get_cube_corner(grid_, _cidx, 0 );
		auto corner1 = TGrid_Traits::get_cube_corner(grid_, _cidx, 1 );
		TGrid_Traits::set_scalar_value( grid_, corner0, 0.5);
		TGrid_Traits::set_scalar_value( grid_, corner1, 0.5);
		
		return;
		for (int i=0; i<8; ++i )
		{
			auto corner = TGrid_Traits::get_cube_corner(grid_, _cidx, i );
			TGrid_Traits::set_scalar_value( grid_, corner, -0.5);
		}
	}
	
	bool starting_cube(unsigned int _cidx)
	{
		auto corner0 = TGrid_Traits::get_cube_corner(grid_, _cidx, 0 );
		auto corner1 = TGrid_Traits::get_cube_corner(grid_, _cidx, 1 );
		auto a = TGrid_Traits::scalar_value( grid_, corner0);
		auto b = TGrid_Traits::scalar_value( grid_, corner1);
		
		return (a == -0.5) && (b == 0.5);
		
		return !TGrid_Traits::is_inside( grid_, corner0 ) && TGrid_Traits::is_inside( grid_, corner1 );
	}
	
	bool ending_cube(unsigned int _cidx)
	{
		auto corner0 = TGrid_Traits::get_cube_corner(grid_, _cidx, 0 );
		auto corner1 = TGrid_Traits::get_cube_corner(grid_, _cidx, 1 );
		auto a = TGrid_Traits::scalar_value( grid_, corner0);
		auto b = TGrid_Traits::scalar_value( grid_, corner1);
		
		return (a == 0.5) && (b == -0.5);
		
		return TGrid_Traits::is_inside( grid_, corner0 ) && !TGrid_Traits::is_inside( grid_, corner1 );
	}
	
	bool boundary_cube(unsigned int _cidx)
	{
		return starting_cube(_cidx) || ending_cube(_cidx);
		Point_descriptor	corner[8];
		unsigned int		cubetype( 0 );
		int i;

		bool opp_bound = false;
		
		for ( i=0; i<8; ++i )
			corner[i] = TGrid_Traits::get_cube_corner(grid_, _cidx, i );

		return TGrid_Traits::is_inside( grid_, corner[0] ) ^ TGrid_Traits::is_inside( grid_, corner[1] );
		
		for (i=0; i<8; ++i )
		{
			if (TGrid_Traits::scalar_value( grid_, corner[i]) == -0.5)
			{
				opp_bound = true;
			}
			if ( !TGrid_Traits::is_inside( grid_, corner[i] ))
			{
				cubetype |= ( 1<<i );
			}
		}

		if ( cubetype == 0 || cubetype == 255 )
			return false;
		else
			return true && opp_bound;
	}
};

}//algorithms

}//hmira

#endif //__VOXELIZE_HPP__