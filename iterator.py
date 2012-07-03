#Iterator.py
#Author: Matias Estrada

from subprocess import Popen
import os
import subprocess



SEED = 12803
ITER = " 1000"
N_ROUTES = " 15"
N_RESTART = " 30"
VECIND = " 30"
N = 30


def main():

    f = open('./sol/output.txt', 'w+')

    for i in xrange(1, N):        
        process = Popen(["./urtp", "td1.txt", "tt1.txt", N_ROUTES, ITER, N_RESTART, VECIND, " " + str(SEED*i) ],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print SEED*i
        
        #print process.stdout
        buff = ""
        for line in process.stdout:
            
            stat, fo, nr = line.split(',')
            if stat == "FINISH":
                f.write( buff + fo + " & "+ nr.rstrip('\n') + r' \\ \hline ' + '\n')
            else:
                buff = fo + " & "+ nr.rstrip('\n') + "& "
        #replace(r"\n", ' ')
        exit_code = os.waitpid(process.pid, 0)
        output = process.communicate()[0]

    f.close()


if __name__ == "__main__":
    main()