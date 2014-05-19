from PIL import Image
import glob
import sys

def to_greyscale(full_filename):
	prefix = full_filename[:-4]
	postfix = full_filename[-4:]
	Image.open(full_filename).convert('L').save(prefix+"_gray"+postfix)

if __name__ == '__main__':
	if len(sys.argv) > 1:
	    filename = sys.argv[1]
	    to_greyscale(filename)
	else:
	    for filename in glob.glob("*.png"):
	            print filename
	            to_greyscale(filename)
