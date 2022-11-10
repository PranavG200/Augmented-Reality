import numpy as np

def est_homography(X, Y):
    """
    Calculates the homography H of two planes such that Y ~ H*X
    If you want to use this function for hw5, you need to figure out
    what X and Y should be.
    Input:
        X: 4x2 matrix of (x,y) coordinates
        Y: 4x2 matrix of (x,y) coordinates
    Returns:
        H: 3x3 homogeneours transformation matrix s.t. Y ~ H*X

    """

    ##### STUDENT CODE START #####
    A =[]
    for i in range(4):
      A.append([-X[i,0], -X[i,1], -1, 0, 0, 0, X[i,0]*Y[i,0], X[i,1]*Y[i,0], Y[i,0]])
      A.append([0, 0, 0, -X[i,0], -X[i,1], -1, X[i,0]*Y[i,1], X[i,1]*Y[i,1], Y[i,1]])

    A = np.array(A)
    [U, S , Vt ] = np.linalg.svd(A)

    H = np.transpose(Vt)[:,-1]
    H = np.reshape(H,(3,3))

    ##### STUDENT CODE END #####

    return H
