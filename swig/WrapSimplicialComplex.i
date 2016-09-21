typedef void CASC;

CASC* newSimplicialComplex();
void deleteSimplicialComplex(CASC* casc);
void insertVertex(CASC* casc, unsigned int n, double x, double y, double z);
void insertFace(CASC* casc, unsigned int a, unsigned int b, unsigned int c);
