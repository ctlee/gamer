#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "SurfaceMesh.h"

void print_vertices(const SurfaceMesh& mesh){
    for(auto x : mesh.get_level<1>()) {
        std::cout << x << ", ";
    }
    std::cout << std::endl;
}

void print_faces(const SurfaceMesh& mesh){
	for(auto x : mesh.get_level_id<3>()) {
		auto name = mesh.get_name(x);
		auto data = mesh.get(x);

		std::cout 	<< "Face(id1=" 	<< name[0]
					<< ", id2=" 	<< name[1]
					<< ", id3="		<< name[2]
					<< ", marker="	<< data.marker
					<< ", selected=" << data.selected
					<< ")";
	}
}


void generateHistogram(const SurfaceMesh& mesh){
    std::array<double,18> histogram;
    histogram.fill(0);

  	for(auto face : mesh.get_level_id<3>()) {
  		auto vertexIDs = mesh.get_name(face);
  		// Unpack the ID's for convenience	
  		Vertex a = mesh.get<1>({vertexIDs[0]}); 
  		Vertex b = mesh.get<1>({vertexIDs[1]}); 
  		Vertex c = mesh.get<1>({vertexIDs[2]});

  		auto binAngle = [&](double angle) -> int{
  			return std::floor(angle/10);
  		};
  		histogram[binAngle(angle(a,b,c))]++;
  		histogram[binAngle(angle(b,a,c))]++;
  		histogram[binAngle(angle(c,a,b))]++;
  		
  	} 
    const int factor = mesh.size<3>()*3;
    std::for_each(histogram.begin(), histogram.end(), [&factor](double& n){
            n = 100.0*n/factor;});

  	for (int x=0; x< 18; x++)
  		std::cout << x*10 << "-" << (x+1)*10 << ": " << std::setprecision(2)  << std::fixed << histogram[x] << std::endl;
  	std::cout << std::endl << std::endl;
}

void translate(SurfaceMesh& mesh, Vector v){
    for(auto& vertex : mesh.get_level<1>())
        vertex += v;
}

void translate(SurfaceMesh& mesh, double dx, double dy, double dz){ 
    Vector v = Vector();
    v[0] = dx; v[1] = dy; v[2] = dz;
    translate(mesh, v);
}

void scale(SurfaceMesh& mesh, Vector v){
    for (auto& vertex : mesh.get_level<1>()){
        vertex[0] *= v[0];
        vertex[1] *= v[1];
        vertex[2] *= v[2];
    }
}

void scale(SurfaceMesh& mesh, double sx, double sy, double sz){
    Vector v = Vector();
    v[0] = sx; v[1] = sy; v[2] = sz;
    scale(mesh, v);
}

void scale(SurfaceMesh& mesh, double s){
    Vector v = Vector();
    v[0] = s; v[1] = s; v[2] = s;
    scale(mesh, v);
}

bool smoothMesh(const SurfaceMesh &mesh, std::size_t minAngle, std::size_t maxAngle, std::size_t maxIter, bool preserveRidges){
	return false;
}