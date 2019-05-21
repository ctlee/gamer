import pytest
import pygamer

class TestFace(object):
    def test_constructors(self):
        face = pygamer.Face(1,2,True)

        assert face.orientation == 1
        assert face.marker == 2
        assert face.selected == True

        face = pygamer.Face()

        assert face.orientation == 0
        assert face.marker == 0
        assert face.selected == False


        face = pygamer.Face(3, True)

        assert face.orientation == 0
        assert face.marker == 3
        assert face.selected == True


class TestSurfaceMesh(object):
    def test_assignment(self):
        mesh = pygamer.SurfaceMesh()
        mesh.insert([1,2,3])

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

        assert pygamer.Face(1,2,True) == fid.data()