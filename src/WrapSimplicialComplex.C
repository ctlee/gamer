#include "SimplicialComplex.h"
#include "Vertex.h"
#include "util.h"
#include "Orientable.h"
#include "WrapSimplicialComplex.h"

struct FaceProperties
{
    int m;
    bool sel;
};

struct Face : Orientable, FaceProperties {
    Face() {}
    Face(Orientable orient, FaceProperties prop) : Orientable(orient), FaceProperties(prop) {}
};

struct Global
{
    bool  closed;                /**< @brief is the surface mesh closed or not */
    int   _marker;               /**< @brief doman marker, to be used when tetrahedralizing */
    float volume_constraint;     /**< @brief volume constraint of the tetrahedralized domain */
    bool  use_volume_constraint; /**< @brief flag that determines if the volume constraint is used */

    float min[3];                /**< @brief minimal coordinate of nodes */
    float max[3];                /**< @brief maximal coordinate of nodes */

    float avglen;                /**< @brief average edge length */

    bool hole;                   /**< @brief flag that determines if the mesh is a hole or not */
};

struct complex_traits
{
    using KeyType = size_t;
    using NodeTypes = util::type_holder<Global,Vertex,void,Face>;
    using EdgeTypes = util::type_holder<Orientable,Orientable,Orientable>;
};

using ASC = simplicial_complex<complex_traits>; 
using WSC = WrappedSimplicialComplex;

WSC::WrappedSimplicialComplex(){
    ASC* asc =  new ASC();
    this->ptr_asc = (void*) asc;
    }
WSC::~WrappedSimplicialComplex(){
    ASC* asc = (ASC*) this->ptr_asc; 
    delete asc;
    }

void WSC::insertVertex(size_t n, double x, double y, double z){
    ASC* asc = (ASC*) this->ptr_asc; 
    Vertex v = Vertex(x,y,z);
    asc->insert<1>({n},v);
}

void WSC::insertFace(size_t a, size_t b, size_t c){
    ASC* asc = (ASC*) this->ptr_asc; 
    asc->insert<3>({a,b,c});
}

void WSC::removeVertex(size_t n){
    ASC* asc = (ASC*) this->ptr_asc; 
    asc->remove<1>({n});
}

void WSC::removeFace(size_t a, size_t b, size_t c){
    ASC* asc = (ASC*) this->ptr_asc;
    asc->remove<3>({a,b,c});
}

int WSC::numVertices() const{
    ASC* asc = (ASC*) this->ptr_asc; 
    return asc->size<1>(); 
}

int WSC::numFaces() const{
    ASC* asc = (ASC*) this->ptr_asc; 
    return asc->size<3>();
}

std::string WSC::as_string() const{
    ASC* asc = (ASC*) this->ptr_asc;
    std::string out = "";
    for(auto& x : asc->get_level<1>()){
        out += x.to_string() + ", ";
    }
    return out;
}

void WSC::print(){
    ASC* asc = (ASC*) this->ptr_asc;
    for(auto& x : asc->get_level<1>()){
        std::cout << x << ", ";
    }
    std::cout << std::endl;
}
