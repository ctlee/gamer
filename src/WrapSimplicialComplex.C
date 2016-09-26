#include "SimplicialComplex.h"
#include "Util.h"
#include "WrapSimplicialComplex.h"

struct Vec3 {
    double x;
    double y;
    double z;
};

struct Orientable {
    int orientation;
};

struct Face : Orientable {
    Face(int orient, long m) : Orientable{orient}, mark(m) {}
    long mark;
};

struct complex_traits
{
    using KeyType = size_t;
    using NodeTypes = util::type_holder<int,Vec3,int, Orientable>;
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

template<> void WSC::insertVertex<>(complex_traits::KeyType n, 
        double x, double y, double z){
    ASC* asc = (ASC*) this->ptr_asc; 
    Vec3 v{x,y,z};
    asc->insert<1>({n},v);
}

template<> void WSC::insertFace<>(complex_traits::KeyType a, 
        complex_traits::KeyType b, complex_traits::KeyType c){
    ASC* asc = (ASC*) this->ptr_asc; 
    asc->insert<3>({a,b,c});
}

int WSC::numVertices(){
    ASC* asc = (ASC*) this->ptr_asc; 
    return asc->size<1>(); 
}

int WSC::numFaces(){
    ASC* asc = (ASC*) this->ptr_asc; 
    return asc->size<3>();
}
