

For Blender 2.79b

conda create --name py35 python=3.5


conda install --file requirements.txt

sudo apt install python3-pip python3-numpy
pip3 install -r requirements.txt
python3 setup.py build -- -DBUILD_BLENDER=on


## Authors
**[Christopher Lee](https://github.com/ctlee)**<br />
Department of Chemistry & Biochemistry<br />
University of California, San Diego

### Contributors to GAMer
* Zeyun Yu (UCSD) and Yuhui Cheng(UCSD)<br />
Development of GAMer v1. To acknowledge your use of GAMer 1, please cite:<br />
[Yu, Z.; Holst, M. J.; Cheng, Y.; McCammon, J. A. Feature-Preserving Adaptive Mesh Generation for Molecular Shape Modeling and Simulation. J. Mol. Graph. Model. 2008, 26 (8), 1370–1380.](https://doi.org/10.1016/j.jmgm.2008.01.007)

* Tom Bartol (Salk Institute) and Johan Hake<br />
Development of Blender GAMer addon.

See also the list of [contributors](https://github.com/ctlee/gamer/contributors) who have contributed to this project.