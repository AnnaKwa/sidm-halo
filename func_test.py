arrX=[0,1,3,7,10]
arrY=[1,1,1,2,3]

def func(x,y):
    z=[]
    for i in range(0,len(x)-1):
        z.append(arrX[i+1]-arrX[i])
    for i in range(0,len(x)-1):
        arrX[i]=arrX[i]+z[i]


func(arrX,arrY)

for i in arrX:
    print i
