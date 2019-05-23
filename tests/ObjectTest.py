import pytest
import pygamer

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
    def test_assignment(self):
        mesh = pygamer.surfacemesh.SurfaceMesh()
        mesh.insertFace([1,2,3])

        fid = mesh.get_simplex_up([1,2,3])
        data = fid.data()

        assert data.orientation == 0
        assert data.marker == 0
        assert data.selected == False

        data.orientation = 1
        data.marker = 2
        data.selected = True

        fid = mesh.get_simplex_up([1,2,3])

        assert data.orientation == 1
        assert data.marker == 2
        assert data.selected == True


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
