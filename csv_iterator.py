#Iterator.py
#Author: Matias Estrada

from subprocess import Popen
import os
import subprocess



SEED = 12803
ITER = " 1000"
N_NODES = " 50"
N_RESTART = " 30"
VECIND = " 30"
N = 5


def main():
    try:

        f = open('./sol/output4.txt', 'w+')

        for i in xrange(16, 20):        
            process = Popen(["./urtp", "td2.txt", "tt2.txt", N_NODES, ITER, N_RESTART, VECIND, " " + str(SEED*i) ],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print SEED*i
            
            #print process.stdout
            buff = ""
            for line in process.stdout:
                
                stat, fo, nr = line.split(',')
                if stat == "FINISH":
                    f.write( buff + fo + ","+ nr.rstrip('\n') + '\n')

            #replace(r"\n", ' ')
            exit_code = os.waitpid(process.pid, 0)
            output = process.communicate()[0]
        f.close()    
    except:
        f.close()


if __name__ == "__main__":
    main()