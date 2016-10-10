#pragma once

#include "biom.h"
#include <tuple>
#include <string>
#include "SurfaceMesh.h"

/** @brief Other data structure SurfaceMeshOld (for surface meshes) */
class SurfaceMeshOld {
public:

    typedef FLTVECT *vertex_iterator;

    // Surface mesh constructors
    SurfaceMeshOld() {}

    SurfaceMeshOld(unsigned int, unsigned int);
    ~SurfaceMeshOld();

    // Move as many of these as possible to global functions.
    // SurfaceMeshOld 
    static SurfaceMeshOld* merge(SurfaceMeshOld *, SurfaceMeshOld *);
    static SurfaceMeshOld* readOFF(std::string filename);
    static SurfaceMeshOld* readPoly(std::string filename);
    static std::tuple<SurfaceMeshOld*,SurfaceMesh_ASC*> readPDB_molsurf(std::string filename);
    static std::tuple<SurfaceMeshOld*,SurfaceMesh_ASC*> readPDB_gauss(const char *filename, float blobbyness, float iso_value);
    static SurfaceMeshOld* sphere(int);
    static std::tuple<SurfaceMeshOld*, SurfaceMesh_ASC*> marchingCube(int, int, int, float *, float, SPNT **);
    static std::tuple<SurfaceMeshOld*, SurfaceMesh_ASC*> marchingCube(int, int, int, float *, float, float *, float, SPNT **);
    static SurfaceMeshOld* readLattice(const char *, float, bool);
    static SurfaceMeshOld* readLattice(const char *, const char *, float, bool);
    static SurfaceMeshOld* triangulate(REAL *pointlist, int numberofpoints, const char *triangle_params);

public:

    const FLTVECT& vertex_attr(int k) const {
        return vertex[k];
    }

    const FLTVECT& vertex_attr(int k, const FLTVECT& v) {
        return vertex[k] = v;
    }

    auto vertices() {
        return util::make_range(vertex, vertex + num_vertices);
    }

    const INT3VECT& face_attr(int k) const {
        return face[k];
    }

    const INT3VECT& face_attr(int k, const INT3VECT& v) {
        return face[k] = v;
    }

    auto faces() {
        return util::make_range(face, face + num_faces);
    }

    int numVertices() const {
        return num_vertices;
    }

    int numFaces() const {
        return num_faces;
    }

    bool isClosed() const {
        return closed;
    }

    bool isHole() const {
        return hole;
    }

    int marker() const {
        return _marker;
    }

    float volumeConstraint() const {
        return volume_constraint;
    }

    bool useVolumeConstraint() const {
        return use_volume_constraint;
    }

public:

    // Almost all of these should *not* be public member functions.
    // We need to re-factor the interface of SurfaceMeshOld so that these can be implemented cleanly
    // in terms of that interface.
    // SurfaceMeshOld should merely handle the proper construction and destruction of a Surface Mesh.
    // But... baby steps... baby steps.

    // Having a needless friend function like this is a special kind of terrible, but one step at a time.
    friend void SurfaceExtract(TeTraMesh *volmesh, SurfaceMeshOld *surfmesh);

    void        createNeighborlist();
    void        destroyNeighborlist();
    void        flipNormals();
    //void        releaseData();
    void        correctNormals();
    void        writeOFF(char *filename);
    void        writePoly(char *filename);
    void        getMinMaxAngles(float *, float *, int *, int *,
                                int, int);
    void        eigenvalues();
    void        refine();
    bool        smooth(unsigned int, unsigned int,
                       unsigned int, bool);
    void        normalSmooth();
    char        coarse(float, float, float, float);
    void        assignActiveSites(ATOM *, unsigned int, unsigned int *);
    ATOM        getCenterRadius(); // this function probably does not belong as an interface function.
    void        translate(float dx, float dy, float dz);
    void        scale(float scale_x, float scale_y, float scale_z);
    void        scale(float scale);
    void        centeralize();
    void        splitMultipleConnectedSurfaces();
    void        removeUnconnectedPatches(int minimal_number);
    void        removeUnconnectedVertices();
    void        deleteVertices();
    void        deleteFaces();

public:

    // these should not be public
    INT3VECT *neighbor;      /**< @brief pointer to the neighbors (triangles) */
    NPNT3   **neighbor_list; /**< @brief pointer to the neighbor list */

private:

    bool  closed;                /**< @brief is the surface mesh closed or not */
    int   _marker;               /**< @brief doman marker, to be used when tetrahedralizing */
    float volume_constraint;     /**< @brief volume constraint of the tetrahedralized domain */
    bool  use_volume_constraint; /**< @brief flag that determines if the volume constraint is used */

    float min[3];                /**< @brief minimal coordinate of nodes */
    float max[3];                /**< @brief maximal coordinate of nodes */

    float avglen;                /**< @brief average edge length */

    int num_vertices;            /**< @brief number of vertices */
    FLTVECT *vertex;             /**< @brief pointer to the vertices */

    bool hole;                   /**< @brief flag that determines if the mesh is a hole or not */

    int num_faces;               /**< @brief number of triangles */
    INT3VECT *face;              /**< @brief pointer to the triangles */
};

