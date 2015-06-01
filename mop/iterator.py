#Iterator.py

__author__ = 'Matias Estrada'
__version__ = '0.1'

from subprocess import Popen
import os
import subprocess



# ./urtp 1 15 5 2 10 0.3 100 10 12345 0.5

# inputid nodes n_routes min_len max_len mut_prob n_iter pob_size seed threshold


SEED = 12803
ITER = (1000, 5000, 10000)
INPUTID = "1"
N_NODES = "15"
MIN_LEN = "2"
MAX_LEN = "7"
N_ROUTES = (3, 4, 5)
MUT_PROB = (0.001, 0.002, 0.005, 0.01, 0.02, 0.05)
POP_SIZE = (5, 10, 20, 50)
THRESHOLD = (0.25, 0.5, 1.0)


def main():

    print "Executing 1800+ instances"

    try:

        f = open('./sol/output_mop.txt', 'w')

        counter = 1

        for niter in ITER:
            for mprob in MUT_PROB:
                for popsize in POP_SIZE:
                    for n_routes in N_ROUTES:
                        for th in THRESHOLD:
                            for i in range(0, 3):
                                print "Execution {0}".format(counter)
                                params = [INPUTID, N_NODES, n_routes, MIN_LEN,
                                MAX_LEN, mprob, niter, popsize, SEED * i, th]
                                params = map(lambda x: str(x), params)
                                process = Popen(["./urtp", ] + params,
                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                                buff = ""
                                for line in process.stdout:
                                        f.write(line)
                                exit_code = os.waitpid(process.pid, 0)
                                output = process.communicate()[0]
                                counter += 1
        f.close()
    except:
        f.close()
        raise


if __name__ == "__main__":
    main()