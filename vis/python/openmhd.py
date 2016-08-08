def data_read(it):
    import numpy as np

    filename = "data/field-%05d.dat" % it
    f = open(filename, 'rb')
    record_marker = np.int32
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=np.double,count=1)
    t0 = buf[0]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=np.int32,count=1)
    ix0 = buf[0]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=np.int32,count=1)
    jx0 = buf[0]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    print( ' t = ', t0 )
    print( ' size = (',ix0,' x ',jx0,')' )

    ix1=0; ix2=ix0-1
    jx1=0; jx2=jx0-1
    
    tmpx = np.zeros((ix0),np.int32)
    tmpy = np.zeros((jx0),np.int32)
    tmp  = np.ndarray((ix0,jx0),np.double)
    data = np.ndarray((ix0,jx0,9),np.double)
    
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmpx = np.fromfile(file=f,dtype=np.double, count=ix0)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmpy = np.fromfile(file=f,dtype=np.double, count=jx0)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = np.fromfile(file=f,dtype=np.double, count=ix0*jx0)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = np.fromfile(file=f,dtype=np.double, count=ix0*jx0)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = np.fromfile(file=f,dtype=np.double, count=ix0*jx0)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = np.fromfile(file=f,dtype=np.double, count=ix0*jx0)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[0:-1,0:-1,4] = tmp[ix1:ix2,jx1:jx2]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[0:-1,0:-1,5] = tmp[ix1:ix2,jx1:jx2]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[0:-1,0:-1,6] = tmp[ix1:ix2,jx1:jx2]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[0:-1,0:-1,7] = tmp[ix1:ix2,jx1:jx2]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[0:-1,0:-1,8] = tmp[ix1:ix2,jx1:jx2]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[0:-1,0:-1,0] = tmp[ix1:ix2,jx1:jx2]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[0:-1,0:-1,1] = tmp[ix1:ix2,jx1:jx2]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[0:-1,0:-1,2] = tmp[ix1:ix2,jx1:jx2]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[0:-1,0:-1,3] = tmp[ix1:ix2,jx1:jx2]
    buf = np.fromfile(file=f,dtype=record_marker,count=1)
    f.close()

    t=t0
    x=tmpx
    y=tmpy
    return x,y,t,data

# end
