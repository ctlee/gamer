import pytest
import pygamer
import math

class TestFace(object):
    def test_constructors(self):
        face = pygamer.surfacemesh.Face(1,2,True)

        assert face.orientation == 1
        assert face.marker == 2
        assert face.selected == True

        face = pygamer.surfacemesh.Face()

        assert face.orientation == 0
        assert face.marker == 0
        assert face.selected == False


        face = pygamer.surfacemesh.Face(3, True)

        assert face.orientation == 0
        assert face.marker == 3
        assert face.selected == True


class TestSurfaceMesh(object):
    def test_insertion(self):
        mesh = pygamer.surfacemesh.SurfaceMesh()
        mesh.insertFace([1,2,3])
        mesh.insertFace([2,3,4])
        assert mesh.nVertices == 4
        assert mesh.nEdges == 5
        assert mesh.nFaces == 2

        mesh.insertEdge([5,6])
        assert mesh.nVertices == 6
        assert mesh.nEdges == 6
        assert mesh.nFaces == 2

        mesh.insertVertex([7])
        assert mesh.nVertices == 7
        assert mesh.nEdges == 6
        assert mesh.nFaces == 2

        mesh.removeVertex([1])
        assert mesh.nVertices == 6
        assert mesh.nEdges == 4
        assert mesh.nFaces == 1


    def test_assignment(self):
        mesh = pygamer.surfacemesh.SurfaceMesh()
        mesh.insertFace([1,2,3])

        fid = mesh.getFace([1,2,3])
        data = fid.data()

        assert data.orientation == 0
        assert data.marker == 0
        assert data.selected == False

        data.orientation = 1
        data.marker = 2
        data.selected = True

        fid = mesh.getFace([1,2,3])
        data = fid.data()

        assert data.orientation == 1
        assert data.marker == 2
        assert data.selected == True

        metadata = mesh.getRoot()
        assert metadata.marker == -1
        assert math.isclose(metadata.volumeConstraint, -1, rel_tol=1e-5) == True
        assert metadata.useVolumeConstraint == False
        assert metadata.ishole == False

        metadata.marker = 5
        metadata.volumeConstraint = 3.14159
        metadata.useVolumeConstraint = True
        metadata.ishole = True

        metadata = mesh.getRoot()
        assert metadata.marker == 5
        assert math.isclose(metadata.volumeConstraint, 3.14159, rel_tol=1e-5) == True
        assert metadata.useVolumeConstraint == True
        assert metadata.ishole == True

        vid = mesh.getVertex([1])
        data = vid.data()
        print(data.position)
        assert data.marker == -1
        assert data.selected == False



class TestVertex(object):
    def test_assignment(self):
        vertex = pygamer.surfacemesh.Vertex()

        vertex[0] = 5
        vertex[1] = 8
        vertex[2] = 10
        vertex.marker = 25
        vertex.selected = True

        assert vertex[0] == 5
        assert vertex[1] == 8
        assert vertex[2] == 10
        assert vertex.marker == 25
        assert vertex.selected == True

        with pytest.raises(IndexError):
            print(vertex[3])
