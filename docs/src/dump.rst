

For Blender 2.79b

conda create --name py35 python=3.5


conda install --file requirements.txt

sudo apt install python3-pip python3-numpy
pip3 install -r requirements.txt
python3 setup.py build -- -DBUILD_BLENDER=on