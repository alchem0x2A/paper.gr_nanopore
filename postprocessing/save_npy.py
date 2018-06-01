import numpy

file_root = ["0001", "0005","001", "005",
             "01", "05", "1" ]

for fr in file_root:
    filename = "./0_{}_out_2Ddata.txt".format(fr)
    print(filename)
    data = numpy.genfromtxt(filename, comments="%")
    npy_name = "0.{}".format(fr)
    numpy.save(npy_name, data)
