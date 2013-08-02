#ifndef MESH_TRAITS_HPP_
#define MESH_TRAITS_HPP_

	template <typename T>
	struct mesh_traits;

	template<typename M>
	struct advanced_mesh_traits;

	struct triangleMesh
	{
		enum { is_triangle = true};
		typedef boost::mpl::true_ is_triangle_t;
	};

	struct polynomialMesh
	{
		enum { is_triangle = false};
		typedef boost::mpl::false_ is_triangle_t;
	};

#endif // MESH_TRAITS_HPP_
