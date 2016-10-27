#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <vector>
#include "SurfaceMesh.h"

void print(const SurfaceMesh & mesh){
    std::cout << "Level: 1" << std::endl;
    for(auto node : mesh.get_level_id<1>()){
        std::cout << "    " << node << std::endl;
    }
    std::cout << "Level: 2" << std::endl;
    for(auto node : mesh.get_level_id<2>()){
        std::cout << "    " << node << std::endl;
    }
    std::cout << "Level: 3" << std::endl;
    for(auto node : mesh.get_level_id<3>()){
        std::cout << "    " << node << std::endl;
    }
}

void print_vertices(const SurfaceMesh& mesh){
    for(auto x : mesh.get_level<1>()) {
        std::cout << x << ", ";
    }
    std::cout << std::endl;
}

void print_faces(const SurfaceMesh& mesh){
	for(auto x : mesh.get_level_id<3>()) {
		auto name = mesh.get_name(x);
		auto data = x.data();

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
  		Vertex a = *mesh.get_node_up<1>({vertexIDs[0]});
  		Vertex b = *mesh.get_node_up<1>({vertexIDs[1]});
  		Vertex c = *mesh.get_node_up<1>({vertexIDs[2]});

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

void edgeFlip(SurfaceMesh& mesh, SurfaceMesh::NodeID<2> edgeID, bool preserveRidges){
    // Assuming that the mesh is manifold
    auto name = mesh.get_name(edgeID);
    std::pair<Vertex, Vertex> shared;    
    shared.first = *mesh.get_node_up<1>({name[0]});
    shared.second = *mesh.get_node_up<1>({name[1]});

    std::pair<Vertex, Vertex> notShared;    
    auto up = mesh.get_cover(edgeID);
    if (up.size() > 2){
        //std::cerr << "This edge participates in more than 2 faces. Returning..." << std::endl;
        return;
    }
    else if (up.size() < 2){
        //std::cerr << "This edge participates in fewer than 2 faces. Returning..." << std::endl;
        return;
    }
    notShared.first  = *mesh.get_node_up<1>({up[0]});
    notShared.second = *mesh.get_node_up<1>({up[1]});

    // Add check to see if notShared.first and second are connected.
    if(mesh.exists<2>({up[0], up[1]})){
        std::cerr << "Found a tetrahedron cannot edge flip." << std::endl;
        return;
    }

    auto getMinAngle = [](const Vertex& a, const Vertex& b, const Vertex& c){
        double minAngle = 999; // dummy for now
        double tmp;
        std::vector<Vertex> triangle = {a,b,c};
        for(int i=0; i < 3; i++){
            std::rotate(triangle.begin(),triangle.begin()+i,triangle.end()) ;
            auto it=triangle.begin();
            tmp = angle(*it, *(it+1), *(it+2));
            if(tmp < minAngle) minAngle = tmp;
        }
        return minAngle;
    };

    // Check if we're on a ridge first
    if(preserveRidges){
        auto a = cross(shared.first-shared.second, shared.first-notShared.first);
        auto b = cross(shared.first-notShared.second, shared.first-shared.second);
        auto val = angle(a,b);
        if (val > 60){
            std::cerr << "Found a ridge, won't flip." << std::endl;
            return;
        }
    }

    // Go through all angle combinations
    double tmp;
    double minAngle = getMinAngle(shared.first, shared.second, notShared.first);
    tmp = getMinAngle(shared.first, shared.second, notShared.second);
    if (tmp < minAngle) minAngle = tmp;

    double minAngleFlip = getMinAngle(notShared.first, notShared.second, shared.first);
    tmp = getMinAngle(notShared.first, notShared.second, shared.second);
    if (tmp < minAngleFlip) minAngleFlip = tmp;

    if (minAngleFlip > minAngle){
        mesh.remove<2>({name[0],name[1]});
        mesh.insert<3>({name[0], up[0], up[1]});
        mesh.insert<3>({name[1], up[0], up[1]});
    }
}

int getValence(SurfaceMesh& mesh, SurfaceMesh::NodeID<1> nodeID){
    std::vector<SurfaceMesh::NodeID<1>> vertices;
    neighbors(mesh, nodeID, std::back_inserter(vertices));
    return vertices.size(); 
}

Vector getNormal(SurfaceMesh& mesh, SurfaceMesh::NodeID<3> faceID){
    auto name = mesh.get_name(faceID);
    auto a = *mesh.get_node_up({name[0]});
    auto b = *mesh.get_node_up({name[1]});
    auto c = *mesh.get_node_up({name[2]});
    return cross(a-b, c-b);
}

Vector getNormal(SurfaceMesh& mesh, SurfaceMesh::NodeID<1> vertexID){
    // find all incident triangles
    auto edges = mesh.up(vertexID);
    auto faces = mesh.up(edges);
    Vector normal; 
    for(auto f : faces) {
        normal += getNormal(mesh, f);
    }
    return normal / (double) faces.size();
}
