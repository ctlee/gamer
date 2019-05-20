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