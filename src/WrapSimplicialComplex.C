#include "SimplicialComplex.h"
#include "Util.h"
#include "WrapSimplicialComplex.h"


extern "C" {
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

    CASC* newSimplicialComplex(){
        ASC* asc = new ASC();
        return (CASC*)asc;
    }

    void deleteSimplicialComplex(CASC* casc){
        ASC* asc = (ASC*) casc; 
        delete asc;
    }

    void insertVertex(CASC* casc, unsigned int n, double x, double y, double z){
        ASC* asc = (ASC*) casc;
        Vec3 v{x,y,z};
        asc->insert<1>({n},v);
    }

    void insertFace(CASC* casc, unsigned int a, unsigned int b, unsigned int c){
        ASC* asc = (ASC*) casc;
        asc->insert<3>({a,b,c});
    }
}
