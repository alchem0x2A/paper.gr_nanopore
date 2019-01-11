import numpy
import os, os.path

# result_path = "../result/concentration/1D"
# file_root = ["0001", "0005","001", "005",
             # "01", "05", "1" ]

result_path = "../result/salt/20/"
file_root = ["NaCl", "LiCl", "MgSO4", "CaCl2", "K2SO4", "K3FeCN", "KCl"]
for fr in file_root:
    filename = os.path.join(result_path,
                            "{}_out_1Ddata.txt".format(fr))
    print(filename)
    try:
        data = numpy.genfromtxt(filename, comments="%")
    except Exception as e:
        print(e)
        continue
    npy_name = os.path.join(result_path, "{}".format(fr))
    numpy.save(npy_name, data)
