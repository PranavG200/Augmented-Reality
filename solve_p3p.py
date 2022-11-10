import numpy as np

def P3P(Pc, Pw, K=np.eye(3)):
    """
    Solve Perspective-3-Point problem, given correspondence and intrinsic

    Input:
        Pc: 4x2 numpy array of pixel coordinate of the April tag corners in (x,y) format
        Pw: 4x3 numpy array of world coordinate of the April tag corners in (x,y,z) format
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: (3,) numpy array describing camera translation in the world (t_wc)

    """
    Pc_3d = np.zeros((4,3))
    Pc_3d[0] = [Pc[0,0]-K[-1,0], Pc[0,1]-K[-1,1], (K[0,0]+K[1,1])/2]
    Pc_3d[1] = [Pc[1,0]-K[-1,0], Pc[1,1]-K[-1,1], (K[0,0]+K[1,1])/2]
    Pc_3d[2] = [Pc[2,0]-K[-1,0], Pc[2,1]-K[-1,1], (K[0,0]+K[1,1])/2]
    Pc_3d[3] = [Pc[3,0]-K[-1,0], Pc[3,1]-K[-1,1], (K[0,0]+K[1,1])/2]

    cpts = (np.linalg.inv(K) @ np.hstack((Pc[1:,:], np.ones((np.shape(Pc[1:,:])[0],1)))).T).T

    c1 = np.dot(cpts[1,:]/np.linalg.norm(cpts[1,:]),cpts[2,:]/np.linalg.norm(cpts[2,:]))
    c2 = np.dot(cpts[0,:]/np.linalg.norm(cpts[0,:]),cpts[2,:]/np.linalg.norm(cpts[2,:]))
    c3 = np.dot(cpts[0,:]/np.linalg.norm(cpts[0,:]),cpts[1,:]/np.linalg.norm(cpts[1,:]))
    print([c1,c2,c3])

    d13 = np.linalg.norm(Pw[1]-Pw[3])
    d23 = np.linalg.norm(Pw[2]-Pw[3])
    d12 = np.linalg.norm(Pw[1]-Pw[2])

    print([d13,d23,d12])

    term_acb = (d23**2 - d12**2)/d13**2
    term_acb_2 = (d23**2 + d12**2)/d13**2

    A4 = (term_acb - 1)**2 - 4*d12**2*(c1)**2/(d13**2)
    A3 = 4*( term_acb*(1-term_acb)*c2 - (1-term_acb_2)*c1*c3 + 2*(d12**2/d13**2)*c1**2*c2)
    A2 = 2*(term_acb**2 - 1 + 2*term_acb**2*c2**2 + 2*((d13**2 - d12**2)/d13**2)*c1**2 - 4*term_acb_2*c1*c2*c3 + 2*((d13**2 - d23**2)/d13**2)*c3**2)
    A1 = 4*(-term_acb*(1+term_acb)*c2 + 2*d23**2*c3**2*c2/d13**2 - (1 - term_acb_2)*c1*c3)
    A0 = (1 + term_acb)**2 - 4*d23**2*c3**2/d13**2

    print([A4,A3,A2,A1,A0])
    v = np.roots([A4, A3, A2, A1, A0])

    R = np.zeros((3,3))
    t = np.zeros(3)
    abc = np.Inf

    for i in range(4):
        if(np.isreal(v[i]) and v[i] > 0):
            u = ( (-1 + term_acb)*v[i]**2 - 2*term_acb*c2*v[i] + 1 + term_acb )/ (2*(c3 - v[i]*c1))
            if(u>0):

                d1 = np.sqrt( d23**2 / ( u**2 + v[i]**2 - 2*u*v[i]*c1))
                d2 = u*d1
                d3 = v[i]*d1
                print([d1,d2,d3])
                p = np.vstack((d1*(cpts[0,:]/np.linalg.norm(cpts[0,:])), d2*(cpts[1,:]/np.linalg.norm(cpts[1,:])), d3*(cpts[2,:]/np.linalg.norm(cpts[2,:]))))
                R_,t_ = Procrustes(Pw[1:,:], p)
                Y = K@ (R_@np.transpose(Pw[0,:]) + t_)
                Y = (Y / Y[-1])[:-1]

                dist = np.linalg.norm(Y - Pc[0,:])
                print(dist)
                if(dist < abc):
                    abc = dist
                    R = R_
                    t = t_

    return np.linalg.inv(R), -np.linalg.inv(R)@t



def Procrustes(X, Y):
    """
    Solve Procrustes: Y = RX + t

    Input:
        X: Nx3 numpy array of N points in camera coordinate (returned by your P3P)
        Y: Nx3 numpy array of N points in world coordinate
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: (3,) numpy array describing camera translation in the world (t_wc)

    """
    ##### STUDENT CODE START #####
    Xc = [0,0,0]
    Yc = [0,0,0]
    size = np.shape(X)[0]

    for i in range(size):
      Xc = Xc + X[i,:]
      Yc = Yc + Y[i,:]

    Xc = Xc/size
    Yc = Yc/size
    X1 = np.transpose(X)
    Y1 = np.transpose(Y)

    for i in range(size):
      X1[:,i] = X1[:,i]- Xc
      Y1[:,i] = Y1[:,i] -Yc

    H = np.matmul(Y1,np.transpose(X1))
    [U,S,Vt] = np.linalg.svd(H)

    S1 = [[1,0,0],[0,1,0],[0,0,np.linalg.det(np.matmul(np.transpose(Vt),np.transpose(U)))]]
    R = np.matmul(U,np.matmul(S1,Vt))

    t = Yc - np.matmul(R,Xc)


    ##### STUDENT CODE END #####

    return R, t
