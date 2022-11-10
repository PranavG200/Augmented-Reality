import numpy as np

def est_pixel_world(pixels, R_wc, t_wc, K):
    """
    Estimate the world coordinates of a point given a set of pixel coordinates.
    The points are assumed to lie on the x-y plane in the world.
    Input:
        pixels: N x 2 coordiantes of pixels
        R_wc: (3, 3) Rotation of camera in world
        t_wc: (3, ) translation from world to camera
        K: 3 x 3 camara intrinsics
    Returns:
        Pw: N x 3 points, the world coordinates of pixels
    """

    ##### STUDENT CODE START #####
    Pc = np.ones((np.shape(pixels)[0],3))
    Pw = np.zeros((np.shape(pixels)[0],3))
    Pc[:,0:2] = pixels
    Rcw = np.transpose(R_wc)
    tcw = -Rcw@(t_wc)

    projmatrix = np.zeros((3,3))
    projmatrix[:,[0,1]] = Rcw[:,[0,1]]
    projmatrix[:,2] = tcw
    for i in range(np.shape(pixels)[0]):
      Pw[i,:] = np.matmul(np.linalg.inv(projmatrix),np.matmul(np.linalg.inv(K),Pc[i,:]))
      Pw[i,0] = Pw[i,0]/Pw[i,2]
      Pw[i,1] = Pw[i,1]/Pw[i,2]
      Pw[i,2] = 0

    ##### STUDENT CODE END #####
    return Pw
