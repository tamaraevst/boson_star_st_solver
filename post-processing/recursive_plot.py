import os
import sys, getopt
import numpy as np
import matplotlib.pyplot as plt

def read_and_write(directory_path):
    for filename in os.listdir(directory_path):
        if filename.startswith("phigrav_"): 
            print(os.path.join(directory_path, filename))

            phigrav = np.loadtxt(os.path.join(directory_path, filename))
            plt.plot(phigrav[:,0], phigrav[:,1], label="{}".format(filename))
            plt.legend(loc="upper right")
            plt.ylabel(r"$\varphi$")
            plt.xlabel(r"$r$")
    plt.show()

def main(argv):
   inputdir = ''
   try:
      opts, args = getopt.getopt(argv,"hi",["idir="])
   except getopt.GetoptError:
      print('test.py -i <input directory>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('test.py -i <input directory>')
         sys.exit()
      elif opt in ("-i", "--idir"):
         inputdir = arg
   
   read_and_write(sys.argv[2])


if __name__ == "__main__":
   main(sys.argv[1])
   print("Done!")