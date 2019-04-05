#-----------------------------------------------------------------------
#    data_read routines for Python 3
#-----------------------------------------------------------------------
#     2016/09/17  S. Zenitani  First version
#     2018/05/04  S. Zenitani  Big-endian variant
#-----------------------------------------------------------------------
# This file contains the following subroutines to load the data.
#
#  * data_read(it,ix1=None,ix2=None,jx1=None,jx2=None)
#  * data_read_from_bigendian(it,ix1=None,ix2=None,jx1=None,jx2=None)
#
# Here we assume a python environment on little-endian computers.
# A little-endian-to-big-endian routine will be straightforwardly obtained.
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#     data_read routine
#-----------------------------------------------------------------------
def data_read(arg1,ix1=None,ix2=None,jx1=None,jx2=None):
    import numpy as np

    if type(arg1) is str:
        filename = arg1
    elif type(arg1) is int:
        filename = "data/field-%05d.dat" % arg1

    f = open(filename, 'rb')
    buf = np.fromfile(file=f,dtype=np.double,count=1)
    t0 = buf[0]
    buf = np.fromfile(file=f,dtype=np.int32,count=1)
    ix0 = buf[0]
    buf = np.fromfile(file=f,dtype=np.int32,count=1)
    jx0 = buf[0]
    print( ' t = ', t0 )
    print( ' size = (',ix0,' x ',jx0,')' )

    if ix1 is None:
        ix1 = 0

    if ix2 is None:
        ix2 = ix0-1

    if jx1 is None:
        jx1 = 0

    if jx2 is None:
        jx2 = jx0-1

    ix = ix2-ix1+1
    jx = jx2-jx1+1

    tmpx = np.ndarray((ix0),np.double)
    tmpy = np.ndarray((jx0),np.double)
    tmp  = np.ndarray((ix0,jx0),np.double)
    data = np.ndarray((ix,jx,9),np.double)

    tmpx = np.fromfile(file=f,dtype=np.double, count=ix0)
    tmpy = np.fromfile(file=f,dtype=np.double, count=jx0)

    # conserved variables (U)
    tmp = np.fromfile(file=f,dtype=np.double, count=ix0*jx0)
    tmp = np.fromfile(file=f,dtype=np.double, count=ix0*jx0)
    tmp = np.fromfile(file=f,dtype=np.double, count=ix0*jx0)
    tmp = np.fromfile(file=f,dtype=np.double, count=ix0*jx0)
    
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,4] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,5] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,6] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,7] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,8] = tmp[ix1:ix2+1,jx1:jx2+1]

    # primitive variables (V)
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,0] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,1] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,2] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype=np.double, count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,3] = tmp[ix1:ix2+1,jx1:jx2+1]
    f.close()

    return tmpx[ix1:ix2+1],tmpy[jx1:jx2+1],t0,data


#-----------------------------------------------------------------------
#     the same data_read code to read big-endian data
#-----------------------------------------------------------------------
def data_read_from_bigendian(arg1,ix1=None,ix2=None,jx1=None,jx2=None):

    if type(arg1) is str:
        filename = arg1
    elif type(arg1) is int:
        filename = "data/field-%05d.dat" % arg1

    f = open(filename, 'rb')
    buf = np.fromfile(file=f,dtype='>d',count=1)
    t0 = buf[0]
    buf = np.fromfile(file=f,dtype='>i',count=1)
    ix0 = buf[0]
    buf = np.fromfile(file=f,dtype='>i',count=1)
    jx0 = buf[0]
    print( ' t = ', t0 )
    print( ' size = (',ix0,' x ',jx0,')' )

    if ix1 is None:
        ix1 = 0

    if ix2 is None:
        ix2 = ix0-1

    if jx1 is None:
        jx1 = 0

    if jx2 is None:
        jx2 = jx0-1

    ix = ix2-ix1+1
    jx = jx2-jx1+1

    tmpx = np.ndarray((ix0),np.double)
    tmpy = np.ndarray((jx0),np.double)
    tmp  = np.ndarray((ix0,jx0),np.double)
    data = np.ndarray((ix,jx,9),np.double)

    tmpx = np.fromfile(file=f,dtype='>d', count=ix0)
    tmpy = np.fromfile(file=f,dtype='>d', count=jx0)

    # conserved variables (U)
    tmp = np.fromfile(file=f,dtype='>d', count=ix0*jx0)
    tmp = np.fromfile(file=f,dtype='>d', count=ix0*jx0)
    tmp = np.fromfile(file=f,dtype='>d', count=ix0*jx0)
    tmp = np.fromfile(file=f,dtype='>d', count=ix0*jx0)

    tmp = (np.fromfile(file=f,dtype='>d', count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,4] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype='>d', count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,5] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype='>d', count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,6] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype='>d', count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,7] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype='>d', count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,8] = tmp[ix1:ix2+1,jx1:jx2+1]

    # primitive variables (V)
    tmp = (np.fromfile(file=f,dtype='>d', count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,0] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype='>d', count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,1] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype='>d', count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,2] = tmp[ix1:ix2+1,jx1:jx2+1]
    tmp = (np.fromfile(file=f,dtype='>d', count=ix0*jx0)).reshape((ix0,jx0),order='F')
    data[:,:,3] = tmp[ix1:ix2+1,jx1:jx2+1]
    f.close()

    return tmpx[ix1:ix2+1],tmpy[jx1:jx2+1],t0,data

# end
#-----------------------------------------------------------------------
