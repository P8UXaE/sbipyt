
import argparse


parser = argparse.ArgumentParser(description="This program does BLA BLA BLA") 
parser.add_argument('-i', '--input',     
                    dest = "infile",                   
                    action = "store",     
                    default = None, 
                    help = "Input file") 
parser.add_argument('-o', '--output',     
                    dest = "outfile",                    
                    action = "store",     
                    default = "default_output.txt",  
                    help = "outputfile") 
options = parser.parse_args()
