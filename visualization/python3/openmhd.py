#-----------------------------------------------------------------------
#    data_read routines for Python 3
#-----------------------------------------------------------------------
#     2016/09/17  S. Zenitani  First version
#     2018/05/04  S. Zenitani  Big-endian variant
#-----------------------------------------------------------------------
# This file contains the following subroutines to load the data.
#
#  * data_read(it,ix1=None,ix2=None,jx1=None,jx2=None,xrange=None,yrange=None)
#  * data_read_from_bigendian(it,ix1=None,ix2=None,jx1=None,jx2=None,xrange=None,yrange=None)
#
# We assume a python environment on little-endian computers.
# A little-endian-to-big-endian version is not provided here, but
# it should be easy to write it.
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#     data_read routine
#-----------------------------------------------------------------------
def data_read(arg1,ix1=None,ix2=None,jx1=None,jx2=None,xrange=None,yrange=None):
    """
    data_read(arg1,ix1=None,ix2=None,jx1=None,jx2=None,xrange=None,yrange=None)

    Reads data from a file.
    The first argument arg1 can be a string or an integer.
    In the string case, the filename can be specified by arg1.
    In the integer case (arg1=N), it reads data from "data/field-0000N.dat".

    Optional arguments
    ------------------
    One can focus on a subdomain by using array indices or tuples.
    ix1 : the first index of a subarray in X (default: 0)
    ix2 : the last index  of a subarray in X (default: nx-1)
    jx1 : the first index of a subarray in Y (default: 0)
    jx2 : the last index  of a subarray in Y (default: ny-1)
    xrange : (x1,x2) indicates the range of x1 < X < x2.
    yrange : (y1,y2) indicates the range of y1 < Y < y2.

    See also
    --------
    To read data from a big-endian file on a little-endian computer,
    use the following routine instead.
    data_read_from_bigendian : It reads data from a big-endian file.
    """
    import numpy as np

    if type(arg1) is int:
        filename = "data/field-%05d.dat" % arg1
    elif type(arg1) is str:
        filename = arg1

    f = open(filename, 'rb')
    buf = np.fromfile(file=f,dtype=np.double,count=1)
    t0 = buf[0]
    buf = np.fromfile(file=f,dtype=np.int32,count=1)
    ix0 = buf[0]
    buf = np.fromfile(file=f,dtype=np.int32,count=1)
    jx0 = buf[0]
    print( ' t = ', t0 )
    print( ' size = (',ix0,' x ',jx0,')' )

    tmpx = np.ndarray((ix0),np.double)
    tmpy = np.ndarray((jx0),np.double)
    tmp  = np.ndarray((ix0,jx0),np.double)
    tmpx = np.fromfile(file=f,dtype=np.double, count=ix0)
    tmpy = np.fromfile(file=f,dtype=np.double, count=jx0)

    if ix1 is None:
        ix1 = 0
    if ix2 is None:
        ix2 = ix0-1
    if jx1 is None:
        jx1 = 0
    if jx2 is None:
        jx2 = jx0-1
    if isinstance(xrange,tuple) or isinstance(xrange,list):
        for i in range(ix0-1,-1,-1):
            if( tmpx[i] < xrange[0] ):
                ix1 = i
                break
        for i in range(0,ix0):
            if( tmpx[i] > xrange[1] ):
                ix2 = i
                break
    if isinstance(yrange,tuple) or isinstance(yrange,list):
        for j in range(jx0-1,-1,-1):
            if( tmpy[j] < yrange[0] ):
                jx1 = j
                break
        for j in range(0,jx0):
            if( tmpy[j] > yrange[1] ):
                jx2 = j
                break

    ix = ix2-ix1+1
    jx = jx2-jx1+1
    data = np.ndarray((ix,jx,9),np.double)

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
#     data_read code to read big-endian data
#-----------------------------------------------------------------------
def data_read_from_bigendian(arg1,ix1=None,ix2=None,jx1=None,jx2=None,xrange=None,yrange=None):
    """
    data_read_from_bigendian(arg1,ix1=None,ix2=None,jx1=None,jx2=None,xrange=None,yrange=None)

    Reads data from a file.
    The first argument arg1 can be a string or an integer.
    In the string case, the filename can be specified by arg1.
    In the integer case (arg1=N), it reads data from "data/field-0000N.dat".

    Optional arguments
    ------------------
    One can focus on a subdomain by using array indices or tuples.
    ix1 : the first index of a subarray in X (default: 0)
    ix2 : the last index  of a subarray in X (default: nx-1)
    jx1 : the first index of a subarray in Y (default: 0)
    jx2 : the last index  of a subarray in Y (default: ny-1)
    xrange : (x1,x2) indicates the range of x1 < X < x2.
    yrange : (y1,y2) indicates the range of y1 < Y < y2.

    See also
    --------
    data_read : It reads data from a file without endian conversion.
    """

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

    tmpx = np.ndarray((ix0),np.double)
    tmpy = np.ndarray((jx0),np.double)
    tmp  = np.ndarray((ix0,jx0),np.double)

    tmpx = np.fromfile(file=f,dtype='>d', count=ix0)
    tmpy = np.fromfile(file=f,dtype='>d', count=jx0)

    if ix1 is None:
        ix1 = 0
    if ix2 is None:
        ix2 = ix0-1
    if jx1 is None:
        jx1 = 0
    if jx2 is None:
        jx2 = jx0-1
    if x1 is not None:
        for i in range(ix0-1,-1,-1):
            if( tmpx[i] < x1 ):
                ix1 = i
                # print( ix1, tmpx[i] )
                break
    if x2 is not None:
        for i in range(0,ix0):
            if( tmpx[i] > x2 ):
                ix2 = i
                # print( ix2, tmpx[i] )
                break
    if y1 is not None:
        for j in range(jx0-1,-1,-1):
            if( tmpy[j] < y1 ):
                jx1 = j
                # print( jx1, tmpy[j] )
                break
    if y2 is not None:
        for j in range(0,jx0):
            if( tmpy[j] > y2 ):
                jx2 = j
                # print( jx2, tmpy[j] )
                break

    ix = ix2-ix1+1
    jx = jx2-jx1+1
    data = np.ndarray((ix,jx,9),np.double)

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
