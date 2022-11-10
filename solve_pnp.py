from est_homography import est_homography
import numpy as np

def PnP(Pc, Pw, K=np.eye(3)):
    """
    Solve Perspective-N-Point problem with collineation assumption, given correspondence and intrinsic

    Input:
        Pc: 4x2 numpy array of pixel coordinate of the April tag corners in (x,y) format
        Pw: 4x3 numpy array of world coordinate of the April tag corners in (x,y,z) format
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: (3, ) numpy array describing camera translation in the world (t_wc)

    """

    ##### STUDENT CODE START #####

    # Homography Approach
    # Following slides: Pose from Projective Transformation
    
    H = est_homography(Pw[:,0:2],Pc)
    H = H/H[-1,-1]
    KinvH = np.matmul(np.linalg.inv(K),H)

    H1 = np.zeros((3,3))
    H1[:,[0,1]] = KinvH[:,[0,1]]
    H1[:,2] = np.cross(KinvH[:,0],KinvH[:,1], axis=0)

    [U,S,Vt] = np.linalg.svd(H1)
    S1 = [[1,0,0],[0,1,0],[0,0,np.linalg.det(np.matmul(U,Vt))]]
    R = np.transpose(np.matmul(U,np.matmul(S1,Vt)))

    t = -R@(KinvH[:,2]/np.linalg.norm(KinvH[:,0]))

    ##### STUDENT CODE END #####

    return R, t
