import sys
from smh import legacy

if __name__=="__main__":
    fname_in  = sys.argv[1]
    fname_out = sys.argv[2]
    legacy.convert_v0_1_to_v0_2(fname_in, fname_out)
