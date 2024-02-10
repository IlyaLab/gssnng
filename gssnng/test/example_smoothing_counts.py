from gssnng import smoothing
import scanpy as sc
import time

def run_smoothing_example():

        print("reading data")
        q = sc.datasets.pbmc3k_processed()

        t0 = time.time()
        print('starting the smOOthing')

        q_list = smoothing.smooth_adata(
                adata=q,
                groupby='louvain',
                smooth_mode='connectivity',
                recompute_neighbors=11,
                cores=8
            )

        t1 = time.time()

        print("Adata List with SMooTHed counts.")
        print("Each is a tuple with groupby category and adata as elements.")
        print(len(q_list))
        for qi in q_list:
                print(qi[1] + "  X size: " + str(qi[0].X.shape))

        print('end time: ' + str(t1))
        print('TOTAL TIME: ' + str(t1-t0))
        print("done")

#if __name__ == '__main__':
run_smoothing_example()
