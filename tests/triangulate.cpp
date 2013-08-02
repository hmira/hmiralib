#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include <hmira/algorithms/marching_cubes.hpp>
#include <hmira/algorithms/fill_holes.hpp>

#include <hmira/meshes/OpenMeshX.hpp>

#include <hmira/grids/ScalarGridT.hh>
#include <hmira/grids/ScalarGrid_traits.h>

#include <hmira/preprocessor/debug_mode.hpp>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

/**
 * \ingroup tests
 * 
 */
int main(int argc, char **argv)
{	
	std::string input_header_filename;
	std::string input_dump_filename;
	std::string output_filename;
	std::string rasterization;
	int x,y,z;
	float x_size, y_size, z_size;
	bool fill;

	po::variables_map vm;
	po::options_description desc("Allowed parameters");
	
	try
	{
		desc.add_options()
		("help,h","produce help message")
		("fill-holes,f","fill holes on the resulting mesh")
		("input-header-file,t", po::value<std::string>(&input_header_filename), "input grid header file")
		("input-dump-file,i", po::value<std::string>(&input_dump_filename), "input grid dump file")
		("output-file,o", po::value<std::string>(&output_filename)->default_value("output.obj"), "output mesh .obj file")
		("rasterize,r", po::value<std::string>(&rasterization)->default_value("full"), "type of rasterization [full|faces|edges]")
		("x-size,X", po::value<float>(&x_size)->default_value(3.f), "X size")
		("y-size,Y", po::value<float>(&y_size)->default_value(3.f), "Y size")
		("z-size,Z", po::value<float>(&z_size)->default_value(3.f), "Z size");
		
		po::positional_options_description p;
		p.add("input-file,i", -1);

		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			std::cerr << desc << std::endl;
			return 0;
		}

	
		if (vm.count("input-header-file"))
		{
		}
		else
		{
			std::cerr << "provide and input file\n" << desc << std::endl;
			return -1;
		}
		
		if (vm.count("input-dump-file"))
		{
		}
		else
		{
			std::cerr << "provide and input file\n" << desc << std::endl;
			return -1;
		}
		fill = (vm.count("fill-holes"));
		
		if (rasterization != "full" && rasterization != "faces" && rasterization != "edges")
		{
			std::cerr << "unknown type of rasterization: \"" << rasterization << "\"\n" << desc << std::endl;
			return -1;
		}
	}
	catch(po::error& e)
	{
		std::cerr << "ERROR: " << e.what() << std::endl;
		return 1;
	}

	ScalarGrid_traits<float, IsoEx::ScalarGridT>::read_header(input_header_filename, x, y, z);
	
	IsoEx::ScalarGridT<float> sg = IsoEx::ScalarGridT<float>(
	OpenMesh::VectorT<float, 3>( 0, 0, 0 ),
	OpenMesh::VectorT<float, 3>( x_size, 0, 0 ),
	OpenMesh::VectorT<float, 3>( 0, y_size, 0 ),
	OpenMesh::VectorT<float, 3>( 0, 0, z_size ),
	x + 1,
	y + 1,
	z + 1);
	
	ScalarGrid_traits<float, IsoEx::ScalarGridT>::read_dump(input_dump_filename, sg);
	std::cerr << "[GRID] : file " << input_dump_filename << " loaded." << std::endl;

	OpenMeshExtended output_mesh;
	auto mc = hmira::algorithms::MarchingCubes<
		IsoEx::ScalarGridT<float>,
		OpenMeshExtended,
		ScalarGrid_traits<float, IsoEx::ScalarGridT>,
		advanced_mesh_traits<OpenMeshExtended>>(sg, output_mesh);
	mc.process();
	
	if (fill)
	{
		hmira::algorithms::fill_holes<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(output_mesh);
	}
	
	if (!OpenMesh::IO::write_mesh(output_mesh, output_filename)) 
	{
		std::cerr << "write error\n";
		exit(1);
	}
	else
	{
		std::cerr << "[MESH] : object: " << output_filename << " written" <<std::endl;
	}


}