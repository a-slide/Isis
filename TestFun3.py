import multiprocessing as mp
import time

def worker(i,q):
    '''stupidly simulates long running process'''
    # task        
    q.put(i)

def listener(filename, q):
    '''listens for messages on the q, writes to file. '''

    f = open(filename, 'wb') 
    while 1:
        m = q.get()
        if m == 'kill':
            f.write('killed')
            break
        f.write(str(m) + '\n')
        f.flush()
    f.close()

def manager(filename):
    #must use Manager queue here, or will not work
    manager = mp.Manager()
    q = manager.Queue()    
    pool = mp.Pool(mp.cpu_count())

    #put listener to work first
    watcher = pool.apply_async(listener, (filename, q, ))

    #fire off workers
    jobs = []
    for i in range(80):
        job = pool.apply_async(worker, (i, q))
        jobs.append(job)

    # collect results from the workers through the pool result queue
    for job in jobs: 
        job.get()

    #now we are done, kill the listener
    q.put('kill')
    pool.close()


def main(filename):
    manager(filename)
