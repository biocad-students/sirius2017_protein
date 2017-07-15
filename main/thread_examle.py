import threading

THREADNUM = 4

def multithread(calcstart,calcstop):
    print("I am thread #",threading.get_ident(),"and my range is ",calcstart,calcstop)

threads = []
length = 1001/THREADNUM
calcpos = 0
for n in range(THREADNUM):
    threads.append(threading.Thread(target=multithread,args = (round(calcpos),round(calcpos+length))))
    calcpos+=length

for thrd in threads:
    thrd.start()
    thrd.join()
