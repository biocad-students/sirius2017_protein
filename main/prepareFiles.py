import os
files = os.listdir('result/')
for x in files:
    try:
        os.rename(x+"0.pdb", "result/"+x+".pdb")
    except:
        pass
