
#pragma once

typedef void CASC;

#ifdef __cplusplus
extern "C"{
#endif

CASC* newSimplicialComplex();
void deleteSimplicialComplex(CASC* casc);
void insertVertex(CASC* casc, unsigned int n, double x, double y, double z);
void insertFace(CASC* casc, unsigned int a, unsigned int b, unsigned int c);
int numVertices(CASC* casc);
int numFaces(CASC* casc);

/*
void insertFaceRGB(CASC* casc, unsigned int a, unsigned int b, unsigned int c, 
        double color_r, double color_g, double color_b);
*/
#ifdef __cplusplus
}
#endif
